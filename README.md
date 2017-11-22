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
