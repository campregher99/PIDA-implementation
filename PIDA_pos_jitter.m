classdef PIDA_pos_jitter < handle
    properties ( Access = private )
        Ts_km1
        up_km1; ui_km1; ud_km1; udd_km1; uff_km1;
        is_sat;
        rf_km1; rf_km2; yf_km1; yf_km2;
        edf_k_m1_; eaf_k_m1_; eaf_k_m2_;
        kp; Ti; Td; Ta; N; M; Tf; b; c; n; saturation;
    end

    methods
        function obj = PIDA_pos_jitter(kp, Ti, Td, Ta, N, M, Tf, b, c, n, saturation)
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
            end
        end

        function initialize(obj, Ts_km1)
            obj.Ts_km1 = Ts_km1;
            obj.ui_km1 = 0;
            obj.is_sat = 0;
            obj.is_sat = 0;
            obj.rf_km1 = 0;
            obj.rf_km2 = 0;
            obj.yf_km1 = 0;
            obj.yf_km2 = 0;
            obj.edf_k_m1_ = 0;
            obj.eaf_k_m1_ = 0;
            obj.eaf_k_m2_ = 0;
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
            obj.ui_km1 = u_k_des - (obj.up_km1 + obj.ud_km1 + obj.udd_km1 + obj.uff_km1);
        end

        function [u_k, up_k, ui_k, ud_k, udd_k] = evaluate(obj, y_k, r_k, Ts_k, uff_k)
            [rf_k, yf_k] = obj.input_filter(r_k, y_k, Ts_k);

            % Proportional
            up_k = obj.kp * (obj.b * rf_k - yf_k);

            % Integral and anti wind-up
            Dui_k = obj.kp * Ts_k / obj.Ti * (rf_k - yf_k);
            if obj.is_sat == 1
                Dui_k = min(Dui_k, 0);
            elseif obj.is_sat == -1
                Dui_k = max(Dui_k, 0);
            end
            ui_k = obj.ui_km1 + Dui_k;

            % Derivative filtering
            [edf_k, edf_k_m1, eddf_k, eddf_k_m1, eddf_k_m2] = obj.derivative_filters(obj.c * rf_k - yf_k, Ts_k);

            % Derivative
            ud_k = obj.kp * obj.Td * (edf_k - edf_k_m1) / Ts_k;

            % Double derivative
            udd_k = 2 * obj.kp * obj.Ta * (obj.Ts_km1 * (eddf_k - eddf_k_m1) + Ts_k * (eddf_k_m2 - eddf_k_m1)) / ...
                (obj.Ts_km1^2 * Ts_k + obj.Ts_km1 * Ts_k^2);

            % Compute control action
            u_k = up_k + ui_k + ud_k + udd_k + uff_k;

            % Apply saturation limits and anti wind-up flag
            if u_k > obj.saturation(2)
                u_k = obj.saturation(2);
                obj.is_sat = 1;
            elseif u_k < obj.saturation(1)
                u_k = obj.saturation(1);
                obj.is_sat = -1;
            else
                obj.is_sat = 0;
            end

            % Update stored values
            obj.up_km1 = up_k;
            obj.ui_km1 = ui_k;
            obj.ud_km1 = ud_k;
            obj.udd_km1 = udd_k;
            obj.uff_km1 = uff_k;
            obj.Ts_km1 = Ts_k;
        end
    end

    methods ( Access = private )
        function [rf_k, yf_k] = input_filter(obj, r_k, y_k, Ts_k)

            % Filtering reference and measurement
            switch obj.n
                case 1
                    rf_k = (Ts_k * r_k + obj.Tf * obj.rf_km1) / (Ts_k + obj.Tf);
                    yf_k = (Ts_k * y_k + obj.Tf * obj.yf_km1) / (Ts_k + obj.Tf);
                case 2
                    rf_k = (r_k * (Ts_k + obj.Ts_km1) * Ts_k * obj.Ts_km1 + ...
                        obj.rf_km1 * ((2 * obj.Ts_km1 * obj.Tf + 2 * obj.Tf^2) * (obj.Ts_km1 + Ts_k)) - ...
                        2 * Ts_k * obj.Tf^2 * obj.rf_km2) / ...
                        (obj.Ts_km1 * ((Ts_k + 2 * obj.Tf) * (obj.Ts_km1 + Ts_k) + 2 * obj.Tf^2));
                    yf_k = (y_k * (Ts_k + obj.Ts_km1) * Ts_k * obj.Ts_km1 + ...
                        obj.yf_km1 * ((2 * obj.Ts_km1 * obj.Tf + 2 * obj.Tf^2) * (obj.Ts_km1 + Ts_k)) - ...
                        2 * Ts_k * obj.Tf^2 * obj.yf_km2) / ...
                        (obj.Ts_km1 * ((Ts_k + 2 * obj.Tf) * (obj.Ts_km1 + Ts_k) + 2 * obj.Tf^2));
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

        function [edf_k, edf_k_m1, eaf_k, eaf_k_m1, eaf_k_m2] = derivative_filters(obj_, ek_, Ts_k_)
            
            % Derivative filter
            edf_k = (Ts_k_ * obj_.N * ek_ + obj_.Td * obj_.edf_k_m1_)/(obj_.Td + Ts_k_ * obj_.N);
            
            % Acceleration filter
            eaf_k = (ek_ * (Ts_k_ + obj_.Ts_km1) * obj_.M^2 * Ts_k_ * obj_.Ts_km1 + ...
                obj_.eaf_k_m1_ * ((2 * obj_.M * obj_.Ts_km1 * obj_.Ta + 2 * obj_.Ta^2) * (obj_.Ts_km1 + Ts_k_)) - ...
                2 * Ts_k_ * obj_.Ta^2 * obj_.eaf_k_m2_)/ ...
                (obj_.Ts_km1 * ((obj_.M^2 * Ts_k_ + 2 * obj_.M * obj_.Ta) * (obj_.Ts_km1 + Ts_k_) + 2 * obj_.Ta^2));
            
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
