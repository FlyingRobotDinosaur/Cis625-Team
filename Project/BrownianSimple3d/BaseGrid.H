//////////////////////////////////////////////////////////////////////
// Grid base class that does just the basics.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BASEGRID_H
#define BASEGRID_H

#include "useful.H"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

class NeighborList {
public:
  double v[3][3][3];
};

class BaseGrid {
  friend class SparseGrid;
 
private:
  // Initialize the variables that get used a lot.
  // Also, allocate the main value array.
  void init() {
    basisInv = basis.inverse();
    nynz = ny*nz;
    size = nx*ny*nz;
    val = new double[size];
  }
public:
  // The most obvious of constructors.
  BaseGrid(Matrix3 basis0, Vector3 origin0, int nx0, int ny0, int nz0) {
    basis = basis0;
    origin = origin0;
    nx = abs(nx0);
    ny = abs(ny0);
    nz = abs(nz0);
    
    init();
    zero();
  }

  // Make an orthogonal grid given the box dimensions and resolution.
  BaseGrid(Vector3 box, double dx) {
    dx = fabs(dx);
    box.x = fabs(box.x);
    box.y = fabs(box.y);
    box.z = fabs(box.z);

    // Tile the grid into the system box.
    // The grid spacing is always a bit smaller than dx.
    nx = int(ceil(box.x/dx));
    ny = int(ceil(box.y/dx));
    nz = int(ceil(box.z/dx));
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;
    basis = Matrix3(box.x/nx, box.y/ny, box.z/nz);
    origin = -0.5*box;

    init();
    zero();
  }

  // Make an orthogonal grid given the box dimensions and resolution.
  BaseGrid(Vector3 box, double dx, double dy, double dz) {
    dx = fabs(dx);
    box.x = fabs(box.x);
    box.y = fabs(box.y);
    box.z = fabs(box.z);

    // Tile the grid into the system box.
    // The grid spacing is always a bit smaller than dx.
    nx = int(ceil(box.x/dx));
    ny = int(ceil(box.y/dy));
    nz = int(ceil(box.z/dz));
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;
    basis = Matrix3(box.x/nx, box.y/ny, box.z/nz);
    origin = -0.5*box;

    init();
    zero();
  }

   // The box gives the system geometry.
  // The grid point numbers define the resolution.
  BaseGrid(Matrix3 box, int nx0, int ny0, int nz0) {
    nx = nx0;
    ny = ny0;
    nz = nz0;

    // Tile the grid into the system box.
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;
    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = -0.5*(box.ex() + box.ey() + box.ez());

    init();
    zero();
  }

  // The box gives the system geometry.
  // dx is the approx. resolution.
  // The grid spacing is always a bit larger than dx.
  BaseGrid(Matrix3 box, Vector3 origin0, double dx) {
    dx = fabs(dx);
    
    // Tile the grid into the system box.
    // The grid spacing is always a bit larger than dx.
    nx = int(floor(box.ex().length()/dx))-1;
    ny = int(floor(box.ey().length()/dx))-1;
    nz = int(floor(box.ez().length()/dx))-1;
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = origin0;

    init();
    zero();
  }

  // The box gives the system geometry.
  // dx is the approx. resolution.
  // The grid spacing is always a bit smaller than dx.
  BaseGrid(Matrix3 box, double dx, double dy, double dz) {
    dx = fabs(dx);
    
    // Tile the grid into the system box.
    // The grid spacing is always a bit smaller than dx.
    nx = int(ceil(box.ex().length()/dx));
    ny = int(ceil(box.ey().length()/dy));
    nz = int(ceil(box.ez().length()/dz));
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
    origin = -0.5*(box.ex() + box.ey() + box.ez());

    init();
    zero();
  }

  // Make an exact copy of a grid.
  BaseGrid(const BaseGrid& g) {
    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    init();
    for (int i = 0; i < size; i++) val[i] = g.val[i];
  }

  BaseGrid& operator=(const BaseGrid& g) {
    delete[] val;

    nx = g.nx;
    ny = g.ny;
    nz = g.nz;
    basis = g.basis;
    origin = g.origin;
    
    init();
    for (int i = 0; i < size; i++) val[i] = g.val[i];

    return *this;
  }

  // Make a copy of a grid, but at a different resolution.
  BaseGrid(const BaseGrid& g, int nx0, int ny0, int nz0) : nx(nx0),  ny(ny0), nz(nz0) {
    if (nx <= 0) nx = 1;
    if (ny <= 0) ny = 1;
    if (nz <= 0) nz = 1;

    // Tile the grid into the box of the template grid.
    Matrix3 box = g.getBox();
    basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);

    origin = g.origin;
    init();

    // Do an interpolation to obtain the values.
    for (int i = 0; i < size; i++) {
      Vector3 r = getPosition(i);
      val[i] = g.interpolatePotential(r);
    }
  }

  // Read a grid from a file.
  BaseGrid(const char* fileName) {
       // Open the file.
    FILE* inp = fopen(fileName,"r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR BaseGrid::BaseGrid Couldn't open file %s.\n",fileName);
      exit(-1);
    }
    //printf("Reading dx file %s...\n", fileName);
    
    size = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    basis = Matrix3(1.0);
    origin = Vector3(0.0);    

    int n = 0;
    double x, y, z;
    char line[STRLEN];
    int p, nRead;
    int deltaCount = 0;
    Vector3 base[3];
    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (len < 2) continue;
      if (line[0] == '#') continue;
    
      if (String::isInt(line[0]) && n < size) {
	// Read grid values.
	nRead = sscanf(line, "%lf %lf %lf", &x, &y, &z);
	if (size > 0) {
	  switch(nRead) {
	  case 1:
	    val[n] = x;
	    n++;
	    if (n != size) {
	      fprintf(stderr,"ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n", fileName);
	      printf("line `%s'\n", line);
	    }
	    break;
	  case 2:
	    val[n] = x;
	    val[n+1] = y;
	    n += 2;
	    if (n != size) {
	      fprintf(stderr,"ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n", fileName);
	      printf("line `%s'\n", line);
	    }
	    break;
	  case 3:
	    val[n] = x;
	    val[n+1] = y;
	    val[n+2] = z;
	    n += 3;
	    break;
	  }
	}
      } else if (len > 5) {
	// Read the grid parameters.
	char start[6];
	for (int i = 0; i < 5; i++) start[i] = line[i];
	start[5] = '\0';

	if(strcmp("origi", start) == 0) {
	  // Get an origin line.
	  p = String::firstSpace(line, STRLEN);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  origin = Vector3(x, y, z);
	  //printf("Origin: %.12g %.12g %.12g\n", x, y, z);
	} else if(strcmp("delta", start) == 0) {
	  // Get a delta matrix line.
	  p = String::firstSpace(line, STRLEN);
	  sscanf(&(line[p+1]), "%lf %lf %lf", &x, &y, &z);
	  base[deltaCount] = Vector3(x, y, z);
	  //printf("Delta %d: %.12g %.12g %.12g\n", deltaCount, x, y, z);
	  if (deltaCount < 2) deltaCount = deltaCount + 1;
	} else if(strcmp("objec", start) == 0) {
	  //printf("%s", line);
	  // Get the system dimensions.
	  if (line[7] != '1') continue;
	  int read = sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &nz);
	  //printf("Size: %d %d %d\n", nx, ny, nz);
	  if (read == 3) {
	    size = nx*ny*nz;
	    nynz = ny*nz;
	    val = new double[size];
	    zero();
	  }
	}
      }
    }
    fclose(inp);

    basis = Matrix3(base[0], base[1], base[2]);
    basisInv = basis.inverse();
    if (size == 0 || n != size) {
      fprintf(stderr,"ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n",fileName);
      printf("declared size: %d, items: %d\n", size, n);
      printf("first value: %10g, final value: %.10g\n", val[0], val[n-1]);
      exit(-1);
    }
  }  
  
  // Write without comments.
  virtual void write(const char* fileName) const {
    write(fileName, "");
  }

  // Writes the grid as a file in the dx format.
  virtual void write(const char* fileName, const char* comments) const {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      fprintf(stderr,"ERROR BaseGrid::write Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    // Write the header.
    fprintf(out, "# %s\n", comments);
    fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "origin %.12g %.12g %.12g\n", origin.x, origin.y, origin.z);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exy, basis.eyy, basis.ezy);
    fprintf(out, "delta %.12g %.12g %.12g\n", basis.exz, basis.eyz, basis.ezz);
    fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
    fprintf(out, "object 3 class array type double rank 0 items %d data follows\n", size);
    
    // Write the data.
    int penultima = 3*(size/3);
    int mod = size - penultima;

    int i;
    for (i = 0; i < penultima; i+=3) {
      fprintf(out, "%.12g %.12g %.12g\n", val[i], val[i+1], val[i+2]);
    }
    if (mod == 1) {
      fprintf(out, "%.12g\n", val[size-1]);
    } else if (mod == 2) {
      fprintf(out, "%.12g %.12g\n", val[size-2], val[size-1]);
    }
    fclose(out);
  }

  // Writes the grid data as a single column in the order:
  // nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
  virtual void writeData(const char* fileName) {
    // Open the file.
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    fprintf(out, "%d\n%d\n%d\n", nx, ny, nz);
    fprintf(out, "%.12g\n%.12g\n%.12g\n", origin.x, origin.y, origin.z);
    fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);
    fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);

    for (int i = 0; i < size; i++) fprintf(out, "%.12g\n", val[i]);
    fclose(out);
  }
 
  // Write the valies in a single column.
  virtual void writePotential(const char* fileName) const {
    FILE* out = fopen(fileName, "w");
    for (int i = 0; i < size; i++) fprintf(out, "%.12g\n", val[i]);
    fclose(out);
  }
  
  virtual ~BaseGrid() {
    delete[] val;
  }

  BaseGrid mult(const BaseGrid& g) {
    for (int i = 0; i < size; i++) val[i] *= g.val[i];
    return *this;
  }

  void zero() {
    for (int i = 0; i < size; i++) val[i] = 0.0;
  }
  
  bool setValue(int j, double v) {
    if (j < 0 || j >= size) return false;
    val[j] = v;
    return true;
  }

  bool setValue(int ix, int iy, int iz, double v) {
    if (ix < 0 || ix >= nx) return false;
    if (iy < 0 || iy >= ny) return false;
    if (iz < 0 || iz >= nz) return false;
    int j = iz + iy*nz + ix*ny*nz;

    val[j] = v;
    return true;
  }

  virtual double getValue(int j) const {
    if (j < 0 || j >= size) return 0.0;
    return val[j];
  }

  virtual double getValue(int ix, int iy, int iz) const {
    if (ix < 0 || ix >= nx) return 0.0;
    if (iy < 0 || iy >= ny) return 0.0;
    if (iz < 0 || iz >= nz) return 0.0;
    
    int j = iz + iy*nz + ix*ny*nz;
    return val[j];
  }

  Vector3 getPosition(int ix, int iy, int iz) const {
    return basis.transform(Vector3(ix, iy, iz)) + origin;
  }

  Vector3 getPosition(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);

    return basis.transform(Vector3(ix, iy, iz)) + origin;
  }

  // Does the point r fall in the grid?
  // Obviously this is without periodic boundary conditions.
  bool inGrid(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);

    if (l.x < 0.0 || l.x >= nx) return false;
    if (l.y < 0.0 || l.y >= ny) return false;
    if (l.z < 0.0 || l.z >= nz) return false;
    return true;
  }

  bool inGridInterp(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);

    if (l.x < 2.0 || l.x >= nx-3.0) return false;
    if (l.y < 2.0 || l.y >= ny-3.0) return false;
    if (l.z < 2.0 || l.z >= nz-3.0) return false;
    return true;
  }

  Vector3 transformTo(Vector3 r) const {
    return basisInv.transform(r-origin);
  }
  Vector3 transformFrom(Vector3 l) const {
    return basis.transform(l) + origin;
  }

  IndexList indices(int j) const {
    int iz = j%nz;
    int iy = (j/nz)%ny;
    int ix = j/(nz*ny);
    IndexList ret;
    ret.add(ix);
    ret.add(iy);
    ret.add(iz);
    return ret;
  }
  int indexX(int j) const { return j/(nz*ny); }
  int indexY(int j) const { return (j/nz)%ny; }
  int indexZ(int j) const { return j%nz; }
  int index(int ix, int iy, int iz) const { return iz + iy*nz + ix*ny*nz; }
  
  int index(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);
    
    int ix = int(floor(l.x));
    int iy = int(floor(l.y));
    int iz = int(floor(l.z));

    ix = wrap(ix, nx);
    iy = wrap(iy, ny);
    iz = wrap(iz, nz);
    
    return iz + iy*nz + ix*ny*nz;
  }

  int nearestIndex(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);
    
    int ix = int(floor(l.x + 0.5));
    int iy = int(floor(l.y + 0.5));
    int iz = int(floor(l.z + 0.5));

    ix = wrap(ix, nx);
    iy = wrap(iy, ny);
    iz = wrap(iz, nz);
    
    return iz + iy*nz + ix*ny*nz;
  }

  int length() const {
    return size;
  }
  Vector3 getOrigin() const {return origin;}
  Matrix3 getBasis() const {return basis;}
  Matrix3 getInverseBasis() const {return basisInv;}
  int getNx() const {return nx;}
  int getNy() const {return ny;}
  int getNz() const {return nz;}
  int getSize() const {return size;}
  void setBasis(const Matrix3& b) {
    basis = b;
    basisInv = basis.inverse();
  }
  void setOrigin(const Vector3& o) {
    origin = o;
  }
  
  // A matrix defining the basis for the entire system.
  Matrix3 getBox() const {
    return Matrix3(nx*basis.ex(), ny*basis.ey(), nz*basis.ez());
  } 
  // The longest diagonal of the system.
  Vector3 getExtent() const {
    return basis.transform(Vector3(nx,ny,nz));
  }
  // The longest diagonal of the system.
  double getDiagonal() const {
    return getExtent().length();
  }
  // The position farthest from the origin.
  Vector3 getDestination() const {
    return basis.transform(Vector3(nx,ny,nz)) + origin;
  }
  // The center of the grid.
  Vector3 getCenter() const {
    return basis.transform(Vector3(0.5*nx,0.5*ny,0.5*nz)) + origin;
  }
  // The volume of a single cell.
  double getCellVolume() const {
    return fabs(basis.det());
  }
  // The volume of the entire system.
  double getVolume() const {
    return getCellVolume()*size;
  }
  Vector3 getCellDiagonal() const {
    return basis.ex() + basis.ey() + basis.ez();
  }
  // Get the area along x and y.
  double getCellArea() const {
    return basis.ex().cross(basis.ey()).length();
  }
  // Get the area along x and y.
  double getArea() const {
    return getCellArea()*nx*ny;
  }


  // Add a fixed value to the grid.
  void shift(double s) {
    for (int i = 0; i < size; i++) val[i] += s;
  }

  // Multiply the grid by a fixed value.
  void scale(double s) {
    for (int i = 0; i < size; i++) val[i] *= s;
  }

  // Get the mean of the entire grid.
  double mean() const {
    double sum = 0.0;
    for (int i = 0; i < size; i++) sum += val[i];
    return sum/size;
  }

  
  // Compute the average profile along an axis.
  // Assumes that the grid axis with index "axis" is aligned with the world axis of index "axis".
  void averageProfile(const char* fileName, int axis) {
    FILE* out = fopen(fileName,"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName);
      exit(-1);
    }

    int dir0 = wrap(axis, 3);
    int dir1 = (axis+1)%3;
    int dir2 = (axis+2)%3;

    int jump[3];
    jump[0] = nynz;
    jump[1] = nz;
    jump[2] = 1;

    int n[3];
    n[0] = nx;
    n[1] = ny;
    n[2] = nz;
   
    for (int i0 = 0; i0 < n[dir0]; i0++) {
      double sum = 0;

      for (int i1 = 0; i1 < n[dir1]; i1++) {
	for (int i2 = 0; i2 < n[dir2]; i2++) {
	  int j = i0*jump[dir0] + i1*jump[dir1] + i2*jump[dir2];
	  sum += val[j];
	}
      }
      
      double v = sum/(n[dir1]*n[dir2]);
      double x = 0.0;
      switch (dir0) {
      case 0:
	x = origin.x + i0*basis.exx;
	break;
      case 1:
	x = origin.y + i0*basis.eyy;
	break;
      case 2:
	x = origin.z + i0*basis.ezz;
	break;
      }
      fprintf(out, "%0.10g %0.10g\n", x, v);
    }

    fclose(out);
  }

  // Get the potential at the closest node.
  virtual double getPotential(Vector3 pos) const {
    // Find the nearest node.
    int j = nearestIndex(pos);

    return val[j];
  }

  virtual double interpolatePotential(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Get the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0;
  }

  // Added by Rogan for times when simpler calculations are required.
  virtual double interpolatePotentialLinearly(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Get the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[2][2][2];
    for (int ix = 0; ix < 2; ix++) {
      for (int iy = 0; iy < 2; iy++) {
	for (int iz = 0; iz < 2; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }

    // Mix along x.
    double g2[2][2];
    for (int iy = 0; iy < 2; iy++) {
      for (int iz = 0; iz < 2; iz++) {

	// p = w[0] * g[0][iy][iz] + (1-w[0]) * g[1][iy][iz]
	// 

	g2[iy][iz] = (1.0-w[0])*g1[0][iy][iz] + w[0]*g1[1][iy][iz];
      }
    }

    // Mix along y.
    double g3[2];
    for (int iz = 0; iz < 2; iz++) {
      g3[iz] = (1.0-w[1])*g2[0][iz] + w[1]*g2[1][iz];
    }

    // DEBUG
    //printf("(0,0,0)=%.1f (0,0,1)=%.1f (0,1,0)=%.1f (0,1,1)=%.1f (1,0,0)=%.1f (1,0,1)=%.1f (1,1,0)=%.1f (1,1,1)=%.1f ",
    //   g1[0][0][0], g1[0][0][1], g1[0][1][0], g1[0][1][1], g1[1][0][0], g1[1][0][1], g1[1][1][0], g1[1][1][1] );
    //printf ("%.2f\n",(1.0-w[2])*g3[0] + w[2]*g3[1]);

    // Mix along z
    return (1.0-w[2])*g3[0] + w[2]*g3[1];
  }


  double interpolateDiffX(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	//g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
	g2[iy][iz] = 3.0*a3*w[0]*w[0] + 2.0*a2*w[0] + a1;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  }

  double interpolateDiffY(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      //g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
      g3[iz] = 3.0*a3*w[1]*w[1] + 2.0*a2*w[1] + a1;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  }

  double interpolateDiffZ(Vector3 pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);
    int homeX = int(floor(l.x));
    int homeY = int(floor(l.y));
    int homeZ = int(floor(l.z));
    
    // Get the array jumps with shifted indices.
    int jump[3];
    jump[0] = nz*ny;
    jump[1] = nz;
    jump[2] = 1;

    // Shift the indices in the home array.
    int home[3];
    home[0] = homeX;
    home[1] = homeY;
    home[2] = homeZ;

    // Shift the indices in the grid dimensions.
    int g[3];
    g[0] = nx;
    g[1] = ny;
    g[2] = nz;

    // Get the interpolation coordinates.
    double w[3];
    w[0] = l.x - homeX;
    w[1] = l.y - homeY;
    w[2] = l.z - homeZ;

    // Find the values at the neighbors.
    double g1[4][4][4];
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
	for (int iz = 0; iz < 4; iz++) {
	  // Wrap around the periodic boundaries. 
	  int jx = ix-1 + home[0];
	  jx = wrap(jx, g[0]);
	  int jy = iy-1 + home[1];
	  jy = wrap(jy, g[1]);
	  int jz = iz-1 + home[2];
	  jz = wrap(jz, g[2]);
	  
	  int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
	  g1[ix][iy][iz] = val[ind];
	}
      }
    }
    double a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    double g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {

	a3 = 0.5*(-g1[0][iy][iz] + 3*g1[1][iy][iz] - 3*g1[2][iy][iz] + g1[3][iy][iz]);
	a2 = 0.5*(2*g1[0][iy][iz] - 5*g1[1][iy][iz] + 4*g1[2][iy][iz] - g1[3][iy][iz]);
	a1 = 0.5*(-g1[0][iy][iz] + g1[2][iy][iz]);
	a0 = g1[1][iy][iz];

	g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    double g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5*(-g2[0][iz] + 3*g2[1][iz] - 3*g2[2][iz] + g2[3][iz]);
      a2 = 0.5*(2*g2[0][iz] - 5*g2[1][iz] + 4*g2[2][iz] - g2[3][iz]);
      a1 = 0.5*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5*(-g3[0] + 3*g3[1] - 3*g3[2] + g3[3]);
    a2 = 0.5*(2*g3[0] - 5*g3[1] + 4*g3[2] - g3[3]);
    a1 = 0.5*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(3.0*a3*w[2]*w[2] + 2.0*a2*w[2] + a1);
  }

  virtual Vector3 interpolateForce(Vector3 pos) const {
    Vector3 f;
    f.x = interpolateDiffX(pos);
    f.y = interpolateDiffY(pos);
    f.z = interpolateDiffZ(pos);
    Vector3 f1 = basisInv.transpose().transform(f);
    return f1;
  }

  static double wrapFloat(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    return x;
  }

  static double wrapDiff(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    if (x >= 0.5*l) x -= l;
    return x;
  }

  static int wrap(int i, int n) {
    if (i < 0) {
      i %= n;
      i += n;
    }
    // The portion above allows i == n, so no else keyword.
    if (i >= n) i %= n;

    return i;
  }

  Vector3 wrap(Vector3 r) const {
    Vector3 l = basisInv.transform(r-origin);
    l.x = wrapFloat(l.x, nx);
    l.y = wrapFloat(l.y, ny);
    l.z = wrapFloat(l.z, nz);
    return basis.transform(l) + origin;
  }

  Vector3 wrapDiff(Vector3 r) const {
    Vector3 l = basisInv.transform(r);
    l.x = wrapDiff(l.x, nx);
    l.y = wrapDiff(l.y, ny);
    l.z = wrapDiff(l.z, nz);
    return basis.transform(l);
  }

  Vector3 wrapDiffNearest(Vector3 r) const {
    Vector3 l = basisInv.transform(r);
    l.x = wrapDiff(l.x, nx);
    l.y = wrapDiff(l.y, ny);
    l.z = wrapDiff(l.z, nz);
	
    double length2 = basis.transform(l).length2();
	
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
	for (int dz = -1; dz <= 1; dz++) {
	  //if (dx == 0 && dy == 0 && dz == 0) continue;
	  Vector3 tmp = Vector3(l.x+dx*nx, l.y+dy*ny, l.z+dz*nz);
	  if (basis.transform(tmp).length2() < length2) {
	    l = tmp;
	    length2 = basis.transform(l).length2();
	  }
	}
      }
    }
	
    return basis.transform(l);
  }


  // Includes the home node.
  // indexBuffer must have a size of at least 27.
  void getNeighbors(int j, int* indexBuffer) const {
    int jx = indexX(j);
    int jy = indexY(j);
    int jz = indexZ(j);
  
    int k = 0;
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
	for (int iz = -1; iz <= 1; iz++) {
	  int ind = wrap(jz+iz,nz) + nz*wrap(jy+iy,ny) + nynz*wrap(jx+ix,nx);
	  indexBuffer[k] = ind;
	  k++;
	}
      }
    }
  }

  // Includes the home node.
  IndexList getNeighborsNoWrap(int j) const {
    int jx = indexX(j);
    int jy = indexY(j);
    int jz = indexZ(j);
  
    IndexList ret;
    int k = 0;
    for (int ix = -1; ix <= 1; ix++) {
      if (jx+ix < 0 || jx+ix >= nx) continue;
      for (int iy = -1; iy <= 1; iy++) {
	if (jy+iy < 0 || jy+iy >= ny) continue;
	for (int iz = -1; iz <= 1; iz++) {
	  if (jz+iz < 0 || jz+iz >= nz) continue;
	  int ind = jz+iz + nz*(jy+iy) + nynz*(jx+ix);
	  ret.add(ind);
	  k++;
	}
      }
    }
    return ret;
  }
  
  // Get the values at the neighbors of a node.
  // Note that homeX, homeY, and homeZ do not need to be wrapped,
  // since we do it here.
  void getNeighborValues(NeighborList* neigh, int homeX, int homeY, int homeZ) const {
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
	for (int iz = -1; iz <= 1; iz++) {
	  int ind = wrap(homeZ+iz,nz) + nz*wrap(homeY+iy,ny) + nynz*wrap(homeX+ix,nx);
	  neigh->v[ix+1][iy+1][iz+1] = val[ind];
	}
      }
    }
  }

 
protected:
  Vector3 origin;
  Matrix3 basis;
  int nx, ny, nz;
  int nynz;
  int size;
  Matrix3 basisInv;
  double* val;
  BaseGrid() {}
};
#endif
