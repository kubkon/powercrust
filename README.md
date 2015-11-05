# powercrust

This repo is a git port of the code found on the Nina Amenta's
website:

[http://web.cs.ucdavis.edu/~amenta/powercrust.html](http://web.cs.ucdavis.edu/~amenta/powercrust.html)

The codebase is distributed under GPL license. See
[LICENSE](LICENSE).

## Introduction

The best documentation for these programs are the research papers
about them:

"The power crust", by Nina Amenta, Sunghee Choi and Ravi Kolluri.
Describes the algorithm and the program.
Check our Web sites at the University of Texas at Austin.

Theory about the algorithm (it's guaranteed to give good models on good
input data) is in:
"The power crust, unions of balls and the medial axis transform", which
will appear in the journal Computational Geometry: Theory and Applications.
Also check our Web sites.

## Programs included in this distribution

The `powercrust` program takes a set of 3D points as input and produces
a polygonal surface, the power crust, and an estimate of the medial
axis as a simplicial complex, the power shape.

The `simplify` program is used in powershape simplification.

The `orient` program is a little utility to get the orientation of
the power crust faces.

The `powercrust` program is built on a slightly modified version of
Ken Clarkson's `hull` program
([http://www.netlib.org/voronoi/hull.html](http://www.netlib.org/voronoi/hull.html))
for convex hull, Delaunay triangulation
and Voronoi diagram, and uses Jonathan Shewchuk's robust predicates for
calculating the positions of Voronoi and power diagram vertices.

Some of the output files are in the Geomview OFF file format, which
can be viewed by Geomview, a free geometry visualization program. It
is an ASCII format easy to convert to whatever else you might want.

## Building the programs

In the terminal, run:

```console
$ make
```

To clean up, run:

```console
$ make clean
```

## Running powercrust

The `powercrust` options are:

```console
-t  is cos (pi-alpha), where alpha is the angle b/w "deep" intersecting balls
    (default = 0.4). When a pole gets labeled as inside or outside, it
    propagates the same label to any neighboring pole whose polar ball deeply
    intersects its own. We have not seen any situations in which we wanted to
    change this default value.

-R  estimate for sampling density constant. 1.0 or larger turns it off.
    (default = 0.6). Used to estimate whether Voronoi cells are "well-shaped",
    for handling noise and for the sharp-corners hack. This is the r in
    r-sampling, that is, the minimum distance to the nearest sample on the
    surface, as a fraction of distance to the medial axis. When r is small the
    sampling is expected to be very dense, the Voronoi cells are expected to be
    really long and skinny and the poles of fat Voronoi cells are thrown away.
    On noise-free inputs with no sharp corners, you can give a value >= 1 for
    -R, and perhaps get a good reconstruction from a sparser input sample.
    Setting -R small might help get good reconstructions from dense but noisy
    samples.

-B  throw away both poles for cells which are not long and skinny. Use the -B
    flag on noise-free inputs from surfaces with sharp corners, and make sure
    the value given with -R is less than one (the smaller the -R value, the
    more it will interpolate data near the corners from the data on adjacent
    smooth surfaces).

-D  no propagation for 1st pole of Voronoi cells which are not long and skinny.
    Use the -D flag when reconstructing surfaces with boundaries, it allows a
    sample on the boundary to have two outside poles.

-w  same as -t, but for trying to label unlabeled poles, the second time around
    (default = 0.3). Once you start fooling around with the -D,-R and/or -B
    options, some poles might fail to be labeled by the regular algorithm. This
    parameter is passed to a second-pass clean-up function which should be a
    little more liberal in propagating labels. Make the -w value smaller when
    you see lots of messages about "unlabeled pole".
    Note: the boundaries between the power cells of labeled and
    unlabeled poles are NOT output as part of the power crust.

-p  This option is used when the input file is actually a set of poles along
    with their labels. It skips the pole computation and labeling steps of the
    algorithm. Use this option for recomputing the power crust after you do
    simplification of the power shape.
```

There are also two options which are passed directly to hull:

```console
-i  the name of the input file.
    Input file is just a list of 3D (x,y,z) point coordinates, in ASCII. Can be
    floating point, BUT they will get rounded to integers. To use little
    numbers, less than one, put a big number in the -m option, below. When the
    -p argument is used, the input should be a list of poles, in the format of
    the output file `inpball`: (x,y,z,radius).

-m  multiplier. The first thing hull does is multiply all the floating point
    numbers in the input by this multiplier and round them into integers. If you
    forget the m option, and your input is all fp numbers less than one, you get
    a one-point output.
```

The output is several files:

```console
pc.off       - powercrust
axis.off     - medial axis (powershape)
axisface.off - only triangles of medial axis (powershape)
inpball      - inside polar balls (x y z r for each ball)
inpole       - just centers of inside poles
outpole      - just centers of outside poles
poleinfo     - Input file for ``simplify''. Poles and powershape faces.
```
Fields in the pole spec are `x,y,z,radius,in_or_out,d`.
`in_or_out` says whether the pole belongs to the inner or
outer power shape (1=inner, 2=outer).
`d` is the greatest distance between any two of the four
input samples whose Voronoi vertex is the pole; it is
used by `simplify` for simplification of the power shape.

Note that `*.off` files are Geomview OFF file format, and as such they can be
viewed with any major 3D mesh processing software such as Meshlab.

Example:

```console
$ ./powercrust -m 1000000 -i hotdogs.pts
```

which is equivalent to:

```console
$ ./powercrust -m 1000000 -i hotdogs.pts -t 0.4 -w 0.3 -R 0.6
```

## Using simplify

The `simplify` program is to be used for powershape
simplification. The powershape which is output by running `powercrust`
(in the file `poleinfo`) can be used as an input to this program.  The
output is a set of poles that are left after simplification. The set
of poles can be input to `powercrust` with the `-p` option to get the
simplified powershape and powercrust.

The arguments to simplify are:

```console
-i <filename> Specifies the input file name.
-o <filename> Specifies the output file name, default set to output.
-n <nt>       The noise threshold.
-r <rt>       The redundancy threshold.
-p 0 or 1     By default, simplify outputs a list of poles and the
              maximum distance between the samples on the the surface of the
              polar balls. This file can be used to plot this information and
              for choosing suitable values for the noise threshold give with -n.
              Using 0 suppresses this output.
-u 0 or 1     By default, simplify removes all unlabeled poles in the
              input. Using -u 0 suppresses this criteria for removal.
```

The list of simplified poles is written in the output file. Simplify
also writes a file `distancelog` which contains a sorted list of the
distance parameters of all poles. This information can be used to pick
a suitable value for the noise threshold given with -n.

Example:

```console
$ ./simplify -i poleinfo -o simp_poles -n 1.0 -r 0.3
```

## Using orient

`orient` is used to set the orientation of the faces in the
powercrust. The output of powercrust has most of the faces in the right
orientation. However, some of them might be wrong because of numerical
errors. `orient` assumes that most of the faces in the input surface
are correctly oriented and reverses the orientation of "incorrect" faces.

The arguments to orient are:

```console
-i <filename> The input file.
-o <filename> The output file.
-r            With this argument the orientation of all the faces is reversed.
```

Example:

```console
$ ./orient -i pc.off -o final.off
```
