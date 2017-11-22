# Coriolis Visualizer
*Really messy source code for the tool I made to make a video on the Coriolis effect*

There are many problems with this source code, including:
* Hard-coded variables that really should be command line arguments: FPS, output directory, frame count output format, etc.
* Animations are stored as commented-out sections of code.
* Too many `#define`s scattered all over the place.
* Redundant calculation of some frames.
* No way to actually output a video. The program outputs  a sequence of frames that must be stitched together in some ohter program.
* Probably some bug in there that causes your USB ports to leak radioactive isotopes or something.

That being said, if you can get something vaguely Earth-like to come out of this program, consider yourself lucky. I was going for the popular pressing-deadlines-write-sloppy-but-fast approach when writing this C.

## Credits

The actual part that I wrote, `main.c`, is licensed by me under the [3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause). See the contents of `main.c` for more details.

`stb_image.h` and `stb_image_write.h` are a part of the freely available [stb libraries](https://github.com/nothings/stb).

`plate.png` is a derivative of the file `world.topo.bathy.200411.3x5400x2700.jpg` widely available at [NASA's Visible Earth](https://visibleearth.nasa.gov/view.php?id=73884) page.
> R. Stöckli, E. Vermote, N. Saleous, R. Simmon and D. Herring (2005).
> The Blue Marble Next Generation - A true color earth dataset including seasonal dynamics from MODIS.
> Published by the NASA Earth Observatory.
> Corresponding author: rѕtοсklі(AT)climate.gsfc.nasa.gov

*A copy of the readme included with the Blue Marble Next Generation file mentioned above is included as `docs/bmng.pdf`. Also, I obscured the email address above a bit to try avoid web spiders from adding it to spam lists. Replace the parentheses and their contents with the corresponding symbol needed in all email addresses. Also do not try to "copy and paste" it.*

`starmap.png` is a derivative of the file `starmap_g4k.jpg` widely available at NASA's [Scientific Visualization Studio](https://svs.gsfc.nasa.gov/3895).
> NASA/Goddard Space Flight Center Scientific Visualization Studio.
> Constellation figures based on those developed for the IAU by Alan MacRobert of Sky and Telescope magazine (Roger Sinnott and Rick Fienberg). 

Do whatever you want with `set.png`. If that's not specific enough, I hereby license `set.png` contained in this repository under [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/).
