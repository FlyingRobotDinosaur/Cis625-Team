///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJECTORYWRITER_H
#define TRAJECTORYWRITER_H

#define PDB_TEMPLATE_LINE "ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      SIM"

#include "cstdio"
#include "useful.H"

class TrajectoryWriter {
public:
  const static int formatTraj = 0;
  const static int formatPdb = 1;

  TrajectoryWriter(const char* filePrefix, const char* formatName, Matrix3 box0, int num0, double timestep0, int outputPeriod0, const String* typeName0) 
    : box(box0), num(num0), timestep(timestep0), outputPeriod(outputPeriod0), typeName(typeName0) {
    pdbTemplate = PDB_TEMPLATE_LINE;
    format = getFormatCode(String(formatName));
    makeUnitCell();

    fileName = filePrefix;
    fileName.add(".");
    fileName.add(getFormatName(format));
  }

  ~TrajectoryWriter() {
  }

private:
  Matrix3 box;
  String fileName;
  int format;
  String pdbTemplate;
  double unitCell[6];
  int num;
  double timestep;
  int outputPeriod;
  const String* typeName;

  void makeUnitCell() {
    double pi = 4.0*atan(1.0);

    unitCell[0] = box.ex().length();
    unitCell[1] = box.ey().length();
    unitCell[2] = box.ez().length();
    
    double bc = box.ey().dot(box.ez());
    double ac = box.ex().dot(box.ez());
    double ab = box.ex().dot(box.ey());

    unitCell[3] = acos(bc/unitCell[1]/unitCell[2])/pi*180.0;
    unitCell[4] = acos(ac/unitCell[0]/unitCell[2])/pi*180.0;
    unitCell[5] = acos(ab/unitCell[0]/unitCell[1])/pi*180.0;
  }

public:
  static int getFormatCode(String format) {
    format.lower();
    if (format == String("traj")) return formatTraj;
    return formatPdb;
  }

  static String getFormatName(int formatCode) {
    switch(formatCode) {
    case formatTraj:
      return String("traj");
      break;
    case formatPdb:
    default:
      return String("pdb");
      break;
    }
  }

  void newFile(const Vector3* pos, const int* type, double t, int n) const {
    switch(format) {
    case formatTraj:
      newTraj(pos, type, t, n);
      break;
    case formatPdb:
    default:
      newPdb(fileName, pos, type);
      break;
    }
  }

  void append(const Vector3* pos, const int* type, double t, int n) const {
    switch(format) {
    case formatTraj:
      appendTraj(pos, type, t, n);
      break;
    case formatPdb:
    default:
      appendPdb(pos, type);
      break;
    }
  }

  void newPdb(const char* outFile, const Vector3* pos, const int* type) const {
    char s[128];

    snprintf(s, 128, "CRYST1   %.3f   %.3f   %.3f  90.00  90.00  90.00 P 1           1\n", box.exx, box.eyy, box.ezz);
    String sysLine(s);

    snprintf(s, 128, "REMARK   frameTime %.10g ns\n", outputPeriod*timestep);
    String remarkLine(s);
    
    String line;

    FILE* out = fopen(outFile, "w");
    fprintf(out, "%s", sysLine.val());
    fprintf(out, "%s", remarkLine.val());

    for (int i = 0; i < num; i++) {
      String name = typeName[type[i]];
      line = makePdbLine(pdbTemplate, i, name, i, name, pos[i], 0.0);
      fprintf(out, "%s",  line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }

  void appendPdb(const Vector3* pos, const int* type) const {
    String line;

    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < num; i++) {
      String name = typeName[type[i]];
      line = makePdbLine(pdbTemplate, i, name, i, name, pos[i], 0.0);
      fprintf(out, "%s", line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }

  void newTraj(const Vector3* pos, const int* type, double t, int n) const {
    FILE* out = fopen(fileName, "w");
    fprintf(out, "TIME %.10g\n", t);
    for (int i = 0; i < n; i++) {
      String name = typeName[type[i]];
      fprintf(out, "%.10g %.10g %.10g %.10g %s\n", t, pos[i].x, pos[i].y, pos[i].z, name.cs());
    }
    fprintf(out, "END\n");
    fclose(out);
  }

  void appendTraj(const Vector3* pos, const int* type, double t, int n) const {
    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < n; i++) {
      String name = typeName[type[i]];
      fprintf(out, "%.10g %.10g %.10g %.10g %s\n", t, pos[i].x, pos[i].y, pos[i].z, name.cs());
    }
    fprintf(out, "END\n");
    fclose(out);
  }


  static String makePdbLine(const String& tempLine, int index, const String& segName, int resId, 
			    const String& name, Vector3 r, double beta) {
    char s[128];

    String record("ATOM  ");
    snprintf(s, 128, "     %5i ", index);
    String si = String(s).range(-6,-1);
    if (name.length() == 4) snprintf(s, 128, "%s   ", name.val());
    else snprintf(s, 128, " %s   ", name.val());
    String nam = String(s).range(0,3);
    String temp0 = tempLine.range(16,21);
  
    snprintf(s, 128,"    %d", resId);
    String res = String(s).range(-4,-1);
    String temp1 = tempLine.range(26,29);
  
    snprintf(s, 128, "       %.3f", r.x);
    String sx = String(s).range(-8,-1);
    snprintf(s, 128, "       %.3f", r.y);
    String sy = String(s).range(-8,-1);
    snprintf(s, 128, "       %.3f", r.z);
    String sz = String(s).range(-8,-1);

    String temp2 = tempLine.range(54,59);
    snprintf(s, 128, "    %.2f", beta);
    String bet = String(s).range(-6,-1);
    String temp3 = tempLine.range(66,71);

    snprintf(s, 128, "%s    ", segName.val());
    String seg = String(s).range(0,3);

    String ret(record);
    ret.add(si);
    ret.add(nam);
    ret.add(temp0);
    ret.add(res);
    ret.add(temp1);
    ret.add(sx);
    ret.add(sy);
    ret.add(sz);
    ret.add(temp2);
    ret.add(bet);
    ret.add(temp3);
    ret.add(seg);
  
    return ret;
  }
};
#endif
