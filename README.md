# PIDA-implementation
This repository present a refernce formulation of the Proportional Integrative Derivatuive Accelerative (PIDA) controller called also Proportional Integrative Double Derivative (PIDD) controller.
The **PIDA control law** can be initially written as:

```math
u(t) = k_p \left( e(t) + \frac{1}{T_i} \int_0^t e(\tau) d\tau + T_d \frac{de(t)}{dt} + T_a \frac{d^2 e(t)}{dt^2} \right)
```

Where:

- `k_p` is the proportional gain  
- `T_i` is the integral time constant  
- `T_d` is the derivative time constant  
- `T_a` is the acceleration (second derivative) time constant

The corresponding **transfer function** is:

```math
C(s) = k_p \left( 1 + \frac{1}{T_i s} + T_d s + T_a s^2 \right)
```

This form, however, is a "textbook" representation and is considered **improper**, since it contains more zeros than poles. In practice, the derivative and acceleration actions must be **filtered** to ensure the controller is implementable and robust.

To address this, **first-order and second-order filters** are used for the derivative and acceleration terms, respectively. The modified controller transfer function becomes:

```math
C(s) = k_p \left( 1 + \frac{1}{T_i s} + \frac{T_d s}{(T_d / N)s + 1} + \frac{T_a s^2}{\left((T_a / M)s + 1\right)^2} \right)
```

Where:

- `N` and `M` are tuning parameters defining the low-pass filter pole locations for the derivative and acceleration terms.

Additionally, to improve performance and reduce sensitivity to high-frequency measurement noise, an **input filter** is added. The transfer function of the input filter is:

```math
f_i(s) = \frac{1}{(T_f s + 1)^n}
```

Where:

- `T_f` is the time constant of the input filter  
- `n` is the order of the filter

The input filter ensures smoother control action by attenuating high-frequency components in the reference and feedback signals.

## Discretization

To implement the PIDA controller in software, the continuous-time control law must be **discretized**.

We start from the time-domain expression:

```math
u(t) = k_p \left( e(t) + \frac{1}{T_i} \int_0^t e(\tau) d\tau + T_d \frac{de(t)}{dt} + T_a \frac{d^2 e(t)}{dt^2} \right)
```

### 1. Integral Term (ZOH approximation)

The integral is approximated using the **zero-order hold (ZOH)** or **Riemann sum**:

```math
\int_0^t e(\tau) d\tau \approx t \cdot e(t)
```

### 2. Derivatives

The first and second derivatives are approximated using **polynomial fitting** (geometrical derivatives):

- First-order polynomial (linear fit) for the first derivative
- Second-order polynomial (quadratic fit) for the second derivative

The derivatives are calculated as:

```math
\pi_1(t) = m t + q \quad \Rightarrow \quad \frac{d\pi_1(t)}{dt} = m
```

```math
\pi_2(t) = a t^2 + b t + c \quad \Rightarrow \quad \frac{d^2\pi_2(t)}{dt^2} = 2a
```

Where (from sampled points \( y_1, y_2, y_3 \) at sampling intervals \( t_1, t_2 \)):

```math
m = \frac{y_3 - y_2}{t_2}
```

```math
a = \frac{-(t_1 y_2 - t_2 y_1 - t_1 y_3 + t_2 y_2)}{t_1 t_2 (t_1 + t_2)}
```

This approach allows handling **non-uniform sampling periods** (i.e., jitter).

---

### 3. Discrete Filter Equations

The **first-order filtered derivative** is based on the Laplace transform equation:

```math
f_d(s) = \frac{e_{d,\text{filt}}(s)}{e_{\text{filt}}(s)} \Rightarrow \frac{T_d}{N} \frac{de_{d,\text{filt}}(t)}{dt} + e_{d,\text{filt}}(t) = e_{\text{filt}}(t)
```

The **second-order filtered acceleration** is:

```math
f_a(s) = \frac{e_{a,\text{filt}}(s)}{e(s)} \Rightarrow \frac{T_a^2}{M^2} \frac{d^2 e_{a,\text{filt}}(t)}{dt^2} + 2 \frac{T_a}{M} \frac{d e_{a,\text{filt}}(t)}{dt} + e_{a,\text{filt}}(t) = e_{\text{filt}}(t)
```

Substituting the discrete derivatives from earlier (with equations for `m` and `a`), we obtain the actual implementations of these filters, which are shown in the pseudo-code listings in the paper (Listing 5).

---

### 4. Input Filter

The input filter is discretized similarly, depending on the order `n` of the filter:

```math
f_i(s) = \frac{1}{(T_f s + 1)^n}
```

For practical implementation, `n` is typically chosen from the set `{0, 1, 2}`. The corresponding code for discrete input filtering is provided in Listing 3 of the paper.

---

### 5. Velocity Form Discretization

For the **velocity form**, the control signal is computed incrementally:

```math
u_k = u_{k-1} + \Delta u_{p,k} + \Delta u_{i,k} + \Delta u_{d,k} + \Delta u_{a,k}
```

Where:
- `u_k` and `u_{k-1}` are the current and previous control values
- `Δ(·) = (·)_k - (·)_{k-1}` represents the delta (change) of each control component

This allows efficient computation and supports features like **anti-windup** and **tracking mode** more naturally.

The full velocity form implementation is detailed in Listing 2 of the paper.

## Controller Features

The proposed reference implementation of the PIDA controller includes several practical features commonly found in industrial PID controllers, adapted to the PIDA framework.

---

### 1. Setpoint Weighting

The controller allows for **setpoint weighting**, which adjusts how the reference signal influences the proportional and derivative actions:

- A parameter `b` is used to weight the setpoint in the proportional term (`b ∈ [0, 1]`).
- A parameter `c` is used for the derivative and acceleration terms (`c ∈ {0, 1}`).

This provides greater flexibility in shaping the closed-loop response and reducing overshoot.

---

### 2. Anti-Windup

To prevent integrator windup when the control signal reaches actuator limits, an **anti-windup mechanism** is included.

- This is implemented using **integrator clamping**.
- The integrator is adjusted based on whether the output is saturated and in which direction.

This prevents the integrator from continuing to accumulate error when the actuator is already saturated, improving recovery performance.

---

### 3. Feedforward

A **feedforward term** can be added to the control signal to compensate for known disturbances or reference changes.

- This term is added directly to the control output and helps reduce lag in the system's response.

Feedforward improves disturbance rejection and tracking without relying solely on feedback.

---

### 4. Tracking Mode / Bump less Switching

The controller supports a **tracking mode** (also called **bump less switching**) to allow smooth transitions between manual and automatic modes.

- When engaged, the controller output tracks a specified value (`u_track`), and the internal states are adjusted accordingly.

This feature ensures there is no sudden jump in the control signal when switching modes, which is important for safety and actuator longevity.

---

### 5. Gain Scheduling

The implementation supports **gain scheduling**, allowing the controller parameters to change dynamically based on operating conditions.

- When done at steady state, this does not introduce any bump.
- During transients:
  - In the **velocity form**, switching is smooth.
  - In the **position form**, a small bump may occur due to the presence of error in the proportional and derivative terms.

This feature enhances performance across varying process conditions.

---

### 6. Jitter Handling

The controller is designed to handle **irregular sampling periods (jitter)**.

- The second derivative, filters, and input processing are adapted to account for changes in sampling time.
- This ensures robustness even when the sampling interval is not constant.

Jitter resilience is critical for digital controllers operating on systems where perfect timing cannot be guaranteed.

## Code organization

The code present in this repository is organized in four main classes, one for each possible scenario.

|                   | **Constant sampling time** | **Irregular sampling time** |
| ----------------- | -------------------------- | --------------------------- |
| **Position form** | PIDA_pos.m                 | PIDA_pos_jitter.m           |
| **Velocity form** | PIDA_vel                   | PIDA_vel_jitter             |

Each class is organized in the same structure in order to guarantee ease of usage. Here the public methods to access the class functionalities are summarized:

-  **Constructor** is the method utilized to create a PIDA controller instance, it requires as input parameter 
	- `k_p` is the proportional gain  
	- `T_i` is the integral time constant  
	- `T_d` is the derivative time constant  
	- `T_a` is the acceleration (second derivative) time constant
	- `N` is the constant defining the relative position of the low-pass filter of the derivative action
	- `M` is the constant defining the relative position of the low-pass filter of the acceleration action
	- `T_f` is the time constant of the input filter
	- `b` and `c` are the set-point weighting factor for the proportional and derivatives action 
	- `n` is the order of the input filter defined in the set `{0,1,2}`
	- `saturation` is a vector defined as `[low_bound up_bound]` 
	- `Ts` in the case of constant sampling time, otherwise it is defined at each evaluation
- **initialize** it initialize the controller with all the inner variable equal to 0, in the case of irregular sampling time, an initial `T_s` is required
- **set_parameters** is utilized to set all the parameters given in the constructor for gain-scheduling purposes
- **get_parameters** gives all the parameters of the current implementation of the cotroller
- **set_control_action** it sets the inner state of the controller to match the desired control action `u_k_des` 
- **evaluate** evaluates the control law based on the measurement `y_k`, the reference `r_k`, and the feed forward action `uff_k`, in case of irregular sampling time also the current sampling period has `Ts_k` to be given