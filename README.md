# Coffee Dynamics

Every single morning without fail, I drink a cup of coffee. Not anything fancy. Just [Nescafe](https://www.nescafe.com/us/products/dark-roast-instant-coffee-7-oz-jar/) instant coffee, dry whole milk, and hot water. Very North African style. I've made this exact cup probably a couple thousand times by now, so I know every single step of it by feel. It's all very intuitve to me. I know the ratio of powder to water, I know when to stop frothing, I know exactly what it's supposed to look like when it's done. 

But something happens during the frothing step that I've never really stopped to think about until recently. When I first put the frother in and turn it on, the flow is very rough and turbulent, it's chaotic. You can literally feel the resistance in your hand. But after a few seconds, something changes and shifts. The mixture smoothens out. The resistance changes and it goes from feeling very violent to feeling almost organized. As if this was all by design. 

This is my attempt to answer that question properly. Not just qualitatively, but mechanistically. The geometry here will be my actual mug. The goal is to model the transition from turbulent to laminar flow in a cup of cofee, and to understand every step of that model. 

<br> 

### My Mug
| Parameter |Symbol | Value |
| -------- | -------- | -------- |
| Inner Radius | R<sub>2 | 0.035 m |
| Outer Radius | -- | 0.0425 m |
| Total height | -- | 0.105 m |
| Liquid depth during frothing | H | 0.04 m |

<img width="300" height="300" alt="mugdimen2" src="https://github.com/user-attachments/assets/2800a922-9245-4265-93f2-65ea8b7e46d2" />


<br>

### My Frother
##### https://www.amazon.com/CIRCLE-JOY-Rechargeable-Cappuccino-Chocolate/dp/B0FT3VTD72?th=1
| Parameter |Symbol | Value |
| -------- | -------- | -------- |
| Disc Radius | R<sub>1 | 0.0105 m |
| Speed setting used | -- | Mid (12,000 RPM) |
| Angular velocity | $\Omega$<sub>1 | 1,257 rad/s |
| Liquid depth during frothing | H | 0.04 m |

<br>

### Misc. Geometry
| Parameter |Symbol | Value |
| -------- | -------- | -------- |
| Gap width | d = R<sub>2 - R<sub>1 | 0.0245 m |
| Outer Radius | $\eta = \frac{R_2}{R_1}$ | 0.30 m |

<br>

### How turbulent is this?
The Reynolds number is the first thing I thought of computing. It is a famous yet relatively simple formula that tells you whether a flow is laminar or turbulent by comparing the inertial forces to the viscous forces. A high Re means inertia dominates (turbulent flow). Low Re means viscosity dominates (laminar flow).

For this particular setup, since our flow is rotational, I used the [Rotational Reynolds Number](https://www.sciencedirect.com/topics/engineering/rotational-reynolds-number) expression defined as:

<div align="center">

$$Re = \frac{\Omega_1 R_1 d}{\nu}$$

</div>

where $\nu$ is the kinematic viscosity of water around 85°C, which is $3.3 \cdot 10^{-7}\ \text{m}^2/\text{s}$. So, plugging everything in:

<div align="center">

$$Re = \frac{1257 \cdot 0.0105 \cdot 0.0245}{3.3 \cdot 10^{-7}} \approx 981,000$$

</div>

Typical pipe flow transistions to turbulence around Re ~ 2,300. This flow is more than 400 times past that. This means that the coffee in my mug every morning is in an extremely violent turbulent state the moment my frother touches it.

### What's actually happening.
Having established that at its worst, my coffee flows with a Reynolds number of nearly a million (981,000), equivalent to the flow an aircraft wing experiences during flight, the question then begs itself. How does this slow down? Where does the smooth, laminar feeling come from, and where does the froth come in?

My first thought was the dissolving of the dry powder into the water. This made the most sense to me intuitively. Undissolved powder particles increase the effective viscosity of the mixture, and as they dissolve the fluid becomes more homogeneous, less resistive. And while this is definitely a factor, it's not the main one. The frother is spinning at 1,257 rad/s. Even a fairly thick, clumpy mixture would be torn to shreds almost instantly at that shear rate. The dissolving happens in the first half second. The smooth feeling comes around 10 seconds later. Something else is sustaining the turbulence after the clumps are gone.

It's explained by two things happening simultaneously. 

First, the bulk fluid spins up. At t=0 the fluid is basically stationary and the frother disc is screaming at 1,257 rad/s. That enormous velocity difference is what drives the turbulence. As you froth, momentum transfers from the disc to the fluid. The whole liquid mass begins to rotate with the frother. The relative velocity between disc and fluid drops, and Re drops with it, because Re depends on that relative velocity.

This part is fairly intuitive. But it still doesn't explain the froth.

The froth doesn't come after laminarization. It comes alongside it, and it actually helps cause it. As air gets entrained into the liquid during frothing, foam forms at the surface. Foam is not a normal fluid. It behaves somewhere between a liquid and a soft solid, and its effective viscosity is orders of magnitude higher than water. Higher viscosity means lower Re. The foam layer actively suppresses turbulent fluctuations from the top down, while spin-up suppresses them from the bottom up. They meet in the middle. The smooth feeling and the froth appearing are not sequential. They are coupled. The froth is part of what's killing the turbulence, not a consequence of it.

To put it simply: in the early violent turbulent phase, air bubbles are being entrained and immediately destroyed by turbulent eddies. The flow is too chaotic to let stable foam structure form. It's only as turbulence weakens that bubbles survive long enough to organize into foam. And there's the feedback loop, weakening turbulence allows foam to form, and forming foam further weakens turbulence.


<p align="center">
  <img width="800" alt="coffeeflow" src="https://github.com/user-attachments/assets/1a65be02-0737-40b8-b75c-573cdc6d034e" />
</p>

### Simplifying Navier-Stokes

The Navier-Stokes equations govern this entire problem. They are a set of equations that describe the motion of viscous fluids. A suite of partial differential equations that mathematically balance acceleration against pressure gradients, viscous forces, and body forces. In their Cartesian form, they are expressed as:

$$\rho\left[\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} + w\frac{\partial u}{\partial z}\right] = -\frac{\partial p}{\partial x} + \mu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}\right) + \rho g_x$$

$$\rho\left[\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} + w\frac{\partial v}{\partial z}\right] = -\frac{\partial p}{\partial y} + \mu\left(\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} + \frac{\partial^2 v}{\partial z^2}\right) + \rho g_y$$

$$\rho\left[\frac{\partial w}{\partial t} + u\frac{\partial w}{\partial x} + v\frac{\partial w}{\partial y} + w\frac{\partial w}{\partial z}\right] = -\frac{\partial p}{\partial z} + \mu\left(\frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} + \frac{\partial^2 w}{\partial z^2}\right) + \rho g_z$$

But for this problem our flow is rotational, so we use the cylindrical form of Navier-Stokes, where the dominant motion is azimuthal. The azimuthal component is expressed as:

$$\rho\left(\frac{\partial u_\theta}{\partial t} + u_r\frac{\partial u_\theta}{\partial r} + \frac{u_\theta}{r}\frac{\partial u_\theta}{\partial \theta} + \frac{u_r u_\theta}{r} + u_z\frac{\partial u_\theta}{\partial z}\right) = -\frac{1}{r}\frac{\partial p}{\partial \theta} + \mu\left(\frac{\partial^2 u_\theta}{\partial r^2} + \frac{1}{r}\frac{\partial u_\theta}{\partial r} - \frac{u_\theta}{r^2} + \frac{1}{r^2}\frac{\partial^2 u_\theta}{\partial \theta^2} + \frac{\partial^2 u_\theta}{\partial z^2}\right) + \rho g_\theta$$

This is a nastier form of NS, but one that can be simplified considerably. Thinking of NS as a form of Newton's Second Law (F = ma) makes this process cleaner. The left side contains the inertial forces, the acceleration of the fluid. The right side contains the forces that drive the flow, the pressure gradient, viscous forces, and gravity.

Before simplifying, some assumptions must be made. I plan to work two cases, one where the flow is steady, and one where it is not. Now obviously the whole point of this project is that the flow transitions from turbulent to laminar over about 10 seconds. If we assume steady state, then by definition the flow is already done transitioning. But to solve this differntial equation by hand, the time term has to go. Before the assumptions, a shoutout to the professor who taught me fluid mechanics, [Dr. Nitesh Nama](https://niteshnama.weebly.com/), a brilliant man.

#### Assumptions

| # | Assumption | Consequence |
|---|---|---|
| 1 | Newtonian fluid | $\mu$ = constant |
| 2 | Incompressible | $\rho$ = constant |
| 3 | Neglect gravity | $g$ = 0 |
| 4 | Axisymmetric | $\partial/\partial\theta$ = 0 |
| 5 | 1D flow | $u_z = u_r = 0$, $\partial/\partial z$ = 0 |
| 6 | No azimuthal pressure gradient | $\partial p/\partial\theta$ = 0 |
| 7 | Steady flow *(Case 1 only)* | $\partial/\partial t$ = 0 |
