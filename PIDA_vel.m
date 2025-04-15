classdef PIDA_vel < handle
    properties ( Access = private )
        Ts;
        u_km1; uff_km1;
        is_sat;
        Drf_km1; Drf_km2; Dyf_km1; Dyf_km2;
        rf_km1; rf_km2; yf_km1; yf_km2;
        edf_k_m1_; eaf_k_m1_; eaf_k_m2_;
        kp; Ti; Td; Ta; N; M; Tf; b; c; n; saturation;
        y_km1; r_km1;
    end

    methods
        function obj = PIDA_vel(kp, Ti, Td, Ta, N, M, Tf, b, c, n, saturation, Ts)
            if nargin < 2
                obj = kp;
            else
                obj.kp = kp;
                obj.Ti = Ti;
                obj.Td = Td;
                obj.Ta = Ta;
                obj.N = N;
                obj.M = M;
                obj.Tf = Tf;
                obj.b = b;
                obj.c = c;
                obj.n = n;
                obj.saturation = saturation;
                obj.Ts = Ts;
            end
        end

        function initialize(obj)
            obj.u_km1 = 0;
            obj.uff_km1 = 0;
            obj.is_sat = 0;
            obj.is_sat = 0;
            obj.Drf_km1 = 0;
            obj.Drf_km2 = 0;
            obj.Dyf_km1 = 0;
            obj.Dyf_km2 = 0;
            obj.rf_km1 = 0; 
            obj.rf_km2 = 0; 
            obj.yf_km1 = 0; 
            obj.yf_km2 = 0;
            obj.edf_k_m1_ = 0;
            obj.eaf_k_m1_ = 0;
            obj.eaf_k_m2_ = 0;
            obj.y_km1 = 0;
            obj.r_km1 = 0;
        end
        
        function set_parameters(obj, kp, Ti, Td, Ta, N, M, Tf, b, c, n, saturation)
            obj.kp = kp;
            obj.Ti = Ti;
            obj.Td = Td;
            obj.Ta = Ta;
            obj.N = N;
            obj.M = M;
            obj.Tf = Tf;
            obj.b = b;
            obj.c = c;
            obj.n = n;
            obj.saturation = saturation;
        end

        function [kp, Ti, Td, Ta, N, M, Tf, b, c, n, saturation] = get_parameters(obj)
            kp = obj.kp;
            Ti = obj.Ti;
            Td = obj.Td;
            Ta = obj.Ta;
            N = obj.N;
            M = obj.M;
            Tf = obj.Tf;
            b = obj.b;
            c = obj.c;
            n = obj.n;
            saturation = obj.saturation;
        end
    
        function set_control_action(obj, u_k_des)
            obj.u_km1 = u_k_des;
        end

        function u_k = evaluate(obj, y_k, r_k, uff_k)
            dr_k = r_k - obj.r_km1;
            dy_k = y_k - obj.y_km1;
            [Drf_k, Dyf_k] = obj.del_input_filter(dr_k, dy_k);
            [rf_k, yf_k] = obj.input_filter(r_k, y_k);

            % Proportional
            Dup_k = obj.kp * (obj.b * Drf_k - Dyf_k);

            % Integral and anti wind-up
            Dui_k = obj.kp * obj.Ts / obj.Ti * (rf_k - yf_k);
            if obj.is_sat == 1
                Dui_k = min(Dui_k, 0);
            elseif obj.is_sat == -1
                Dui_k = max(Dui_k, 0);
            end

            % Derivative filtering
            [Dedf_k, Dedf_k_m1, Deddf_k, Deddf_k_m1, Deddf_k_m2] = obj.derivative_filters(obj.c * Drf_k - Dyf_k);

            % Derivative
            Dud_k = obj.kp * obj.Td * (Dedf_k - Dedf_k_m1) / obj.Ts;

            % Double derivative
            Dudd_k = 2 * obj.kp * obj.Ta * (obj.Ts * (Deddf_k - Deddf_k_m1) + obj.Ts * (Deddf_k_m2 - Deddf_k_m1)) / ...
                (2 * obj.Ts^3);
            
            % Feed forward
            Duff_k = uff_k - obj.uff_km1;

            % Compute control action
            u_k_ = obj.u_km1 + Dup_k + Dui_k + Dud_k + Dudd_k + Duff_k;

            % Apply saturation limits and anti wind-up flag
            if u_k_ > obj.saturation(2)
                u_k = obj.saturation(2);
                obj.is_sat = 1;
            elseif u_k_ < obj.saturation(1)
                u_k = obj.saturation(1);
                obj.is_sat = -1;
            else
                u_k = u_k_;
                obj.is_sat = 0;
            end

            % Update stored values
            obj.u_km1 = u_k_;
            obj.uff_km1 = uff_k;
            obj.r_km1 = r_k;
            obj.y_km1 = y_k;
        end
    end

    methods ( Access = private )
        function [Drf_k, Dyf_k] = del_input_filter(obj, Dr_k, Dy_k)

            % Filtering deviation reference and measurement
            switch obj.n
                case 1
                    Drf_k = (obj.Ts * Dr_k + obj.Tf * obj.Drf_km1) / (obj.Ts + obj.Tf);
                    Dyf_k = (obj.Ts * Dy_k + obj.Tf * obj.Dyf_km1) / (obj.Ts + obj.Tf);
                case 2
                    Drf_k = (Dr_k * 2 * obj.Ts^3 + ...
                        obj.Drf_km1 * ((2 * obj.Ts * obj.Tf + 2 * obj.Tf^2) * 2 *obj.Ts) - ...
                        2 * obj.Ts * obj.Tf^2 * obj.Drf_km2) / ...
                        (obj.Ts * ((obj.Ts + 2 * obj.Tf) * 2 * obj.Ts + 2 * obj.Tf^2));
                    Dyf_k = (Dy_k * 2 * obj.Ts^3 + ...
                        obj.Dyf_km1 * ((2 * obj.Ts * obj.Tf + 2 * obj.Tf^2) * 2 * obj.Ts) - ...
                        2 * obj.Ts * obj.Tf^2 * obj.Dyf_km2) / ...
                        (obj.Ts * ((obj.Ts + 2 * obj.Tf) * 2 * obj.Ts + 2 * obj.Tf^2));
                otherwise
                    Drf_k = Dr_k;
                    Dyf_k = Dy_k;
            end

            % Upgrade stored values
            obj.Drf_km2 = obj.Drf_km1;
            obj.Drf_km1 = Drf_k;
            obj.Dyf_km2 = obj.Dyf_km1;
            obj.Dyf_km1 = Dyf_k;
        end

        function [rf_k, yf_k] = input_filter(obj, r_k, y_k)
            
            % Filtering reference and measurement
            switch obj.n
                case 1
                    rf_k = (obj.Ts * r_k + obj.Tf * obj.rf_km1) / (obj.Ts + obj.Tf);
                    yf_k = (obj.Ts * y_k + obj.Tf * obj.yf_km1) / (obj.Ts + obj.Tf);
                case 2
                    rf_k = (r_k * 2 * obj.Ts^3 + ...
                        obj.rf_km1 * ((2 * obj.Ts * obj.Tf + 2 * obj.Tf^2) * 2 *obj.Ts) - ...
                        2 * obj.Ts * obj.Tf^2 * obj.rf_km2) / ...
                        (obj.Ts * ((obj.Ts + 2 * obj.Tf) * 2 * obj.Ts + 2 * obj.Tf^2));
                    yf_k = (y_k * 2 * obj.Ts^3 + ...
                        obj.yf_km1 * ((2 * obj.Ts * obj.Tf + 2 * obj.Tf^2) * 2 * obj.Ts) - ...
                        2 * obj.Ts * obj.Tf^2 * obj.yf_km2) / ...
                        (obj.Ts * ((obj.Ts + 2 * obj.Tf) * 2 * obj.Ts + 2 * obj.Tf^2));
                otherwise
                    rf_k = r_k;
                    yf_k = y_k;
            end

            % Upgrade stored values
            obj.rf_km2 = obj.rf_km1;
            obj.rf_km1 = rf_k;
            obj.yf_km2 = obj.yf_km1;
            obj.yf_km1 = yf_k;
        end

        function [edf_k, edf_k_m1, eaf_k, eaf_k_m1, eaf_k_m2] = derivative_filters(obj_, ek_)
            
            % Derivative filter
            edf_k = (obj_.Ts * obj_.N * ek_ + obj_.Td * obj_.edf_k_m1_)/(obj_.Td + obj_.Ts * obj_.N);

            % Acceleration filter
            eaf_k = (ek_ * 2 * obj_.Ts^3 * obj_.M^2 + ...
                obj_.eaf_k_m1_ * ((2 * obj_.M * obj_.Ts * obj_.Ta + 2 * obj_.Ta^2) * 2 * obj_.Ts) - ...
                2 * obj_.Ts * obj_.Ta^2 * obj_.eaf_k_m2_)/ ...
                (obj_.Ts * ((obj_.M^2 * obj_.Ts + 2 * obj_.M * obj_.Ta) * 2 * obj_.Ts + 2 * obj_.Ta^2));

            % Passing out the variables
            edf_k_m1 = obj_.edf_k_m1_;
            eaf_k_m1 = obj_.eaf_k_m1_;
            eaf_k_m2 = obj_.eaf_k_m2_;

            % Upgrade stored values
            obj_.edf_k_m1_ = edf_k;
            obj_.eaf_k_m2_ = obj_.eaf_k_m1_;
            obj_.eaf_k_m1_ = eaf_k;

        end
    end
end
