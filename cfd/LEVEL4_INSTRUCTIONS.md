# Level 4: filming the real coffee and comparing it with the simulation

The goal is a top-down slow-motion video of the coffee surface while the frother runs, turned into a measured velocity field with PIVlab and compared with the simulation at the same moments.

## What you need

- Your phone, ideally with a 240 fps slow-motion mode.
- A tripod or anything that holds the phone steady, lens pointing straight down, about 25-30 cm above the mug.
- Tracer particles that float: ground black pepper works well on a light milk-coffee surface. Use a pinch, spread thin.
- A ruler or coin in the frame at the same height as the coffee surface, for scale.
- Even light from the side, not from directly above (that causes glare on the surface).

## Filming

1. Make the coffee exactly as usual, but stop before frothing. Let the surface go still.
2. Sprinkle the pepper so it covers the surface evenly.
3. Start recording in slow motion, wait 1 second, then start the frother at mid speed (12,000 RPM) at the usual depth and position (centered, head about 1 cm off the bottom).
4. Keep filming for about 10 seconds.
5. Note roughly when the surface first dips toward the frother and when the froth starts to form.
6. Film two or three takes. The first second after switching on matters most.

## PIVlab (free, MATLAB app or standalone)

1. Export the slow-motion clip as an image sequence (PIVlab can also read the video directly).
2. Load the frames, draw a circular region of interest on the inside of the mug rim, and mask the frother shaft.
3. Calibrate with the ruler or coin (pixels to mm) and set the real time between frames: 1/240 s for 240 fps slow motion, not the playback rate.
4. Use multi-pass FFT with 64 then 32 then 16 pixel windows.
5. Export the vector field for each frame as CSV (x, y, u, v).

## Comparing

```powershell
# simulated surface velocity, same CSV layout as PIVlab, 2 mm below the still surface
python cfd/scripts/extract_surface_velocity.py --depth 0.002

# side-by-side video; --start is when the frother switches on in your clip (seconds)
python cfd/scripts/video_compare.py --real my_coffee.mp4 --start 1.25 --duration 8
```

Things to compare:

- **Swirl speed vs radius**, averaged around the circle, at t = 0.5, 1, 2 and 5 s: the most direct test of the simulated spin-up.
- **When the central dip reaches the frother** (Level 3 predicts this) vs when you see air being drawn down.
- **Whether the surface rotates like a solid body** (speed proportional to radius) or faster near the center.

## Caveats

- The simulation has no foam. Once froth covers the surface, the pepper and the comparison stop being meaningful.
- The model's frother head is an estimate (a solid torus). If your real whisk differs, update `cfd/geometry/params.py` and rerun.
