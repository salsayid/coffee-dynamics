# Coffee Dynamics

Every single morning without fail, I drink a cup of coffee. Not anyhting fancy. Just [Nescafe](https://www.nescafe.com/us/products/dark-roast-instant-coffee-7-oz-jar/) instant coffee, dry whole milk, and hot water. Very North African style. I've made this exact cup probably a couple thousand times by now, so I know every single step of it by feel. It's all very intuitve to me. I konw the ratio of powder to water, I konw when to stop frothing, I know exactly what it's supposed to look like when it's done. 

But something happens during the frothing step that I've never really stopped to think about until recently. When I first put the frother in and turn it on, the flow is very rough and turbulent, it's chaotic. You can literally feel the resistance in your hand. But after a few seconds, something changes and shifts. The mixutre smoothens out. The resistance changes and it goes from feeling very violent to feeling almost organized. As if this was done by design. 

This is my attempt to answer that question properly. Not just qualitatively, but mechanistically. The geometry here will be my actual mug. The goal is to model the transition from turbulent to laminar flow in a cup of cofee, nad to understand every stop of that model. 

<br> 

### My Mug
| Parameter |Symbol | Value |
| -------- | -------- | -------- |
| Inner Radius | R<sub>2 | 0.035 m |
| Outer Radius | -- | 0.0425 m |
| Total height | -- | 0.105 m |
| Liquid depth during frothing | H | 0.04 m |

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

Pipe flow transistions to turbulence around Re ~ 2,300. This flow is more than 400 times past that. This means that the coffee in my mug every morning is in an extremely violent turbulent state the moment my forther touches it.

### What's actually happening.
Having established that at its worst, my coffee flows with a Reynolds number of nearly a million (981,000), equivalent to the flow an aircraft wing experiences during flight, the question then begs itself. How does this slow down? Where does the smooth, laminar feeling come from, and where does the froth come in?

My first thought was the dissolving of the dry powder into the water. This made the most sense to me intuitively. Undissolved powder particles increase the effective viscosity of the mixture, and as they dissolve the fluid becomes more homogeneous, less resistive. And while this is definitely a factor, it's not the main one. The frother is spinning at 1,257 rad/s. Even a fairly thick, clumpy mixture would be torn to shreds almost instantly at that shear rate. The dissolving happens in the first half second. The smooth feeling comes around 10 seconds later. Something else is sustaining the turbulence after the clumps are gone.

It's explained by two things happening simultaneously. 

First, the bulk fluid spins up. At t=0 the fluid is basically stationary and the frother disc is screaming at 1,257 rad/s. That enormous velocity difference is what drives the turbulence. As you froth, momentum transfers from the disc to the fluid. The whole liquid mass begins to rotate with the frother. The relative velocity between disc and fluid drops, and Re drops with it, because Re depends on that relative velocity.

This part is fairly intuitive. But it still doesn't explain the froth.

The froth doesn't come after laminarization. It comes alongside it, and it actually helps cause it. As air gets entrained into the liquid during frothing, foam forms at the surface. Foam is not a normal fluid. It behaves somewhere between a liquid and a soft solid, and its effective viscosity is orders of magnitude higher than water. Higher viscosity means lower Re. The foam layer actively suppresses turbulent fluctuations from the top down, while spin-up suppresses them from the bottom up. They meet in the middle. The smooth feeling and the froth appearing are not sequential. They are coupled. The froth is part of what's killing the turbulence, not a consequence of it.

To put it simply: in the early violent turbulent phase, air bubbles are being entrained and immediately destroyed by turbulent eddies. The flow is too chaotic to let stable foam structure form. It's only as turbulence weakens that bubbles survive long enough to organize into foam. And there's the feedback loop, weakening turbulence allows foam to form, and forming foam further weakens turbulence.
