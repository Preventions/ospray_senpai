#ifndef __CUBEUTILS_HXX__
#define __CUBEUTILS_HXX__ 1

#include <stdio.h>
#include <string>
#include <vector>
#include <limits>
#include "VectorT.hxx"
#include "Matrix44.hxx"
#include "Array3.hxx"
#include "Clock.hxx"
#include "TransferFunction.hxx"
#include <algorithm>

#if 0
#include <omp.h>
#endif


#ifndef IDXCAST
#define IDXCAST(x) (float&)(x)
#endif    


const char atom_str[256][4] = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", 
  "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
  "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
  "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
  "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", 
                               "Uus", "Uuo", "_A", "_B", "_C", "_D", "_E", "_F", "_G", "_H", "_I", "_J", "_K", "E"};


#define cerr_dbg(x)
//#define cerr_dbg(x) cerr << x

#define cerr_dbg2(x) cerr << x

inline vec4 transform_atom(Matrix44& m, vec4& atom)
{
  vec4 dst;
  Matrix44::transform3(dst, m, atom);
  dst[3] = atom[3];
  return dst;
}

struct CubeUtilsSettings
{
  float ball_radius;
  float stick_distance;
  
  float vdw_radius[256];
  float covalent_radius[256];
  
  vec4f atom_colors[256];

  float stick_distance_map[256][256];
  
  //find wavefunction amplitude psi, instead of rho
  bool use_sqrt;
  bool display_uchar_volume;
  bool display_analytical_volume;
  bool build_sticks;
  bool use_fixed_stick_distance;
  bool use_custom_atoms;
	
  int ball_radius_fixed;
  float rdf_clamp;
  
  int max_atoms;
  int abort_atomic_number;
  
  CubeUtilsSettings(CubeUtilsSettings* defaults)
  {
    set_defaults(defaults);
  }
  
  CubeUtilsSettings()
  {
    use_sqrt = false;
    display_uchar_volume = false;
    display_analytical_volume = false;
    build_sticks = true;
    use_fixed_stick_distance = false;
    use_custom_atoms = false;
    
    ball_radius = .25f;
    stick_distance = 1.35f;
    rdf_clamp = 1.5f;
    max_atoms = 0;
    
    abort_atomic_number = 4096;
    
    for(int i=0; i<256; i++)
    {
      atom_colors[i] = vec4f(1,1,1,1);
      vdw_radius[i] = 0.f;
      covalent_radius[i] = 0.f;
    }
    
    //AARONBAD -- add to this as we go
    atom_colors[1] = vec4f(1,1,1, .35);                   //hydrogen
    atom_colors[3] = vec4f(1,0,1, 1.28);			      //lithium
    atom_colors[5] = vec4f(1.0, 0, 0, .66);           //oxygen copy (AARONBAD)
    atom_colors[6] = vec4f(0.4, 0.4, 0.4, .73);           //carbon
    atom_colors[8] = vec4f(1.0, 0, 0, .66);               //oxygen
    atom_colors[13] = vec4f(.5, .5, .5, 1.21);						//aluminum
    atom_colors[14] = vec4f(0.3, 0.3, 1.0, 1.11);         //silicon
    atom_colors[42] = vec4f(.251, .8784, .8157, 1.54);		//molybdenum
    atom_colors[78] = vec4f(1., .65, 0, 1.36);						//platinum

    //custom atoms for Ken-ichi's new Al2O3 system
    atom_colors[120] =  vec4f(.5, .5, .5, 1.21);          //_A (Al)
    atom_colors[121] =  vec4f(1.0, 0, 0, .66);            //_B (O)
    atom_colors[122] =  vec4f(1.0, 0, 0, .66);            //_C (Al)
    atom_colors[123] =  vec4f(.5, .5, .55, 1.21);          //_D (Al)
    atom_colors[124] =  vec4f(1.0, 0, .05, .66);            //_E (O)
    atom_colors[125] =  vec4f(1.0, 0, .05, .66);            //_F (Al)
    atom_colors[126] =  vec4f(.5, .55, .5, 1.21);          //_G (Al)
    atom_colors[127] =  vec4f(1.0, .05, 0, .66);            //_H (O)
    atom_colors[128] =  vec4f(1.0, .05, 0, .66);            //_I (Al)

    atom_colors[255] = vec4(0.f, 0.f, 1.f, 4.f);          //"unknown" or "special" atom

    covalent_radius[3] = 1.28;
    covalent_radius[6] = .73f;     //77 sp3, 73 sp2, 69 sp     //let the user change this if desired
    covalent_radius[8] = .64f;
    covalent_radius[13] = 1.37f;
    covalent_radius[14] = 1.11f;
    covalent_radius[42] = 1.39f;
    covalent_radius[78] = 1.36f;
    
    vdw_radius[1] = 1.2f;
    vdw_radius[3] = 1.82f;
    vdw_radius[6] = 1.7f;
    vdw_radius[7] = 1.55f;
    vdw_radius[8] = 1.52f;
    vdw_radius[9] = 1.47f;
    vdw_radius[13] = 1.84f;
    vdw_radius[14] = 2.1f;
    vdw_radius[42] = 1.75f;
    vdw_radius[78] = 1.75f;

    for(int i=0; i<256; i++)
      for(int j=0; j<256; j++)
	stick_distance_map[i][j] = 0.f;

    set_stick_distance("H", "H", 1);
    set_stick_distance("O", "O", 1.5);
    set_stick_distance("H", "O", 1.3);
    set_stick_distance("Al", "Al", 3.0);
    set_stick_distance("Li", "Li", 3.0);
    set_stick_distance("H", "Al", 2.0);
    set_stick_distance("H", "Li", 2.0);
    set_stick_distance("O", "Al", 2.3);
    set_stick_distance("O", "Li", 2.3);
    set_stick_distance("Al", "Li", 3.3);    
    
    update_atom_radii();
  }
  
  ~CubeUtilsSettings(){}
  
  void set_defaults(CubeUtilsSettings* dcus)
  {
    memcpy(this, dcus, sizeof(CubeUtilsSettings));
  }

  void set_stick_distance(const char* a1, const char* a2, float d)
  {
    int i1 = get_atomic_number(string(a1));
    int i2 = get_atomic_number(string(a2));
    stick_distance_map[i1][i2] = d;
    stick_distance_map[i2][i1] = d;
  }
  
  void update_atom_radii()
  {    
    for(int i=0; i<256; i++)
      //atom_colors[i][3] = sqrt(vdw_radius[i]) * ball_radius;
      atom_colors[i][3] = sqrt(vdw_radius[i]);
  }

  void parse_args(int argc, char* argv[])
  {
    bool update = false;
    for(int i=1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        if (!strcmp(argv[i], "-stick_distance")) 
        {
          stick_distance = atof(argv[++i]);
          update = true;
        }
        if (!strcmp(argv[i], "-use_fixed_stick_distance")) 
        {
          use_fixed_stick_distance = bool(atoi(argv[++i]));
          update = true;
        }
        else if (!strcmp(argv[i], "-ball_radius")) 
        {
          ball_radius = atof(argv[++i]);
          update = true;
        }   
        else if (!strcmp(argv[i], "-rdf_clamp")) 
        {
          rdf_clamp = atof(argv[++i]);
        }               
        else if (!strcmp(argv[i], "-use_sqrt")) 
        {
          use_sqrt = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-build_sticks")) 
        {
          build_sticks = atoi(argv[++i]);
        }        
        else if (!strcmp(argv[i], "-set_vdw_radius")) 
        {
          int atom = get_atomic_number(string(argv[++i]));
          vdw_radius[atom] = atof(argv[++i]);
          update = true;
        }
        else if (!strcmp(argv[i], "-set_covalent_radius")) 
        {
          int atom = get_atomic_number(string(argv[++i]));
          covalent_radius[atom] = atof(argv[++i]);
          update = true;          
        }
        else if (!strcmp(argv[i], "-set_stick_distance")) 
        {
          int atom1 = get_atomic_number(string(argv[++i]));
          int atom2 = get_atomic_number(string(argv[++i]));
          stick_distance_map[atom1][atom2] = atof(argv[++i]);
          update = true;          
        }
        else if (!strcmp(argv[i], "-max_atoms"))
        {
          max_atoms = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-abort_atomic_number"))
        {
          abort_atomic_number = atoi(argv[++i]);
        }    
        else if (!strcmp(argv[i], "-use_custom_atoms"))
        {
          //use the _A, _B etc custom atoms, i.e. add underscores to the .xyz file
          use_custom_atoms = 1;
        }
        else if (!strcmp(argv[i], "-set_atom_colors"))
        {
          int atom = get_atomic_number(string(argv[++i]));
          atom_colors[atom][0] = atof(argv[++i]);
          atom_colors[atom][1] = atof(argv[++i]);
          atom_colors[atom][2] = atof(argv[++i]);
          atom_colors[atom][3] = atof(argv[++i]);     
          update = true;               
        }
        else if (!strcmp(argv[i], "-display_uchar_volume"))
        {
          display_uchar_volume = true;
          display_analytical_volume = false;				
        }
        else if (!strcmp(argv[i], "-display_analytical_volume"))
        {
          display_analytical_volume = true;
          display_uchar_volume = false;
        }
      }
    }
    
    if (update)
      update_atom_radii();
  }
    
  int get_atomic_number(string str)
  {
    for(int i=0; i<119; i++)
    {
      if (strcmp(atom_str[i], str.c_str()) == 0)
        return i;
    }
    return 255;
  }

  float get_covalent_radius(int a)
  {
    if (covalent_radius[a] == 0.f)
    {
      cerr_dbg(  "Error, no specified covalent radius for " << atom_str[a] << ", please put it in CubeUtilsSettings!" << endl );
      covalent_radius[a] = 1.3001f;
    }
    return covalent_radius[a];
  }
  
  
  float get_vdw_radius(int a)
  {
    if (vdw_radius[a] == 0.f)
    {
      cerr_dbg(  "Error, no specified Van der Waals radius for " << atom_str[a] << ", please put it in CubeUtilsSettings!" << endl  ); 
      vdw_radius[a] = 2.0001f;
    }
    return vdw_radius[a];
  }

  float get_stick_distance_squared(int a, int b)
  {
    if (use_fixed_stick_distance)
      return stick_distance * stick_distance;
    
    const float d = get_stick_distance(a,b);
    cerr_dbg("stick_distance(" << a << ", " << b << ") = " << d << endl);
    return d*d;
  }

  float get_stick_distance_squared()
  {
    //const float d = get_stick_distance(a,b);
    const float d = stick_distance;
    return d * d;
  }
  
  
  float get_stick_distance(int a, int b)
  {
    if (stick_distance_map[a][b])
      return stick_distance_map[a][b];
    else
      return ((get_vdw_radius(a) + get_vdw_radius(b)) * .5f * stick_distance);   
    //return (max(get_covalent_radius(a), get_covalent_radius(b)) + max(get_vdw_radius(a), get_vdw_radius(b))) * .5;
  }
  
  float get_bond_cutoff_length(int a, int b)
  {
    return ((get_covalent_radius(a) + get_covalent_radius(b)) * .5f);    
  }
};


using namespace std;

struct BallsMC
{
  BallsMC(){ balls.reserve(2); }
  ~BallsMC(){}
  vector<int> balls;
};

struct BuildMacrocell
{
  BuildMacrocell()
  {
    balls_no_radius.reserve(2);
    balls.reserve(4);
    sticks.reserve(8);
  }
  ~BuildMacrocell(){}
  
  vector<int> balls_no_radius;
  vector<int> balls;
  vector<int> sticks;
};

/*
struct Macrocell
{
  unsigned char min;
  unsigned char max;
  unsigned char pad[2];
  int ball_start;
  int stick_start;
  int stick_end;
};
 */

#ifdef BAS_TEXTURE_BUFFER
typedef int indextype;
#else
typedef float indextype;
#endif


struct CubeData
{
  CubeUtilsSettings* cus;
                
  int voxels_per_macrocell;
  float voxels_per_angstrom;
  int volume_dims_log2;
  int macrocell_dims_log2;
  int xyz_subtype;
	
	bool write_gridded_particles;
	
  float rdf_pad_percent;
  
  int num_indices;
  indextype* indices;
      
  string filename;
  string filebase;
  string name;
  string description;
  string atoms_string;
  string forces_filename;
  
  string nhdr_options_string;

  vector<vec2i> atom_type_counts;
  vector<int>* bond_adjacencies;
  
  //spaces:
  //cell space [0,1]^3 for whatever the unit cell is
  //grid space [0,volume_dims]^3
  //world space [X, Y] in Angstroms

  //cell space to grid space (must be affine)
  Matrix44 cs2gs_transform;
  Matrix44 gs2cs_transform;

  //world space to grid space (can be non-affine)
  Matrix44 ws2gs_transform;
  Matrix44 gs2ws_transform;
  
  //vectors for .xyz/.raw files where ws2gs is affine. 
  //we use these for generating the ws2gs_transform
  vec4 ws2gs_translate;
  vec4 ws2gs_scale;

  //cell space to world space (can be non-affine)
  Matrix44 ws2cs_transform;
  Matrix44 cs2ws_transform;
  
  vector<vec4f> cs_atoms;     //used by chgcar
  vector<vec4f> gs_atoms;     //on the GPU
  vector<vec4f> ws_atoms;     //used by cube, xyz
  
  vector<vec4f> forces;

  //DEFAULT: we always want to have ws_atoms, and ws2gs_transform
  
  vector<int> special_atoms;
  
  vec3i volume_dims;
  vector<Array3Base*> volumes;
  int cvi;
  
  //Array3<float> float_volume;
  //Array3<unsigned char> uchar_volume;
  
  Array3<vec4f> macrocells;   //pertains to ball & stick
  
  vec3i mc_dims;
 
  vec3f bounds_min;
  vec3f bounds_max;
  
  struct RadialDistanceFunction
  {
    CubeData* cd;
    CubeUtilsSettings* cus;
    int atomid;
    int atom;

    float Z;
    float clamp_radius;
    float sigma;
    float inv_sigma_squared;
    
    virtual void init(CubeData* _cd, CubeUtilsSettings* _cus, int _atomid, int _atom)
    {
      cd = _cd;
      cus = _cus;
      atom = _atom;
      atomid = _atomid;
      Z = sqrtf(float(atom));
      clamp_radius = cus->get_vdw_radius(atom) * cus->rdf_clamp;
      sigma = cus->get_covalent_radius(atom);
      inv_sigma_squared = 1.f / (sigma * sigma);
    }
    
    virtual float get_clamp_radius()
    {
      return clamp_radius;
    }
    
    virtual float evaluate(float distance_squared)
    {
      //note, the 2.506 number is sqrt(2*pi)
      const float r2 = distance_squared * inv_sigma_squared;
      
      return (Z * expf(-r2));

      //for testing with Attila
      //return .125f*Z*exp(-125.f*r2) + Z*exp(-2.f*r2);

      //default
      //return (Z / ( sigma * sigma_mul * 2.5066282746f )) * exp(-1.5f * r2));
    }
  };

  #if 0
  struct RDF_FromDistribution : public RadialDistanceFunction
  {
    RadialDensityDistribution* rdd;
    
    virtual float evaluate(float distance)
    {
      return density_data[c][i]
    
      //note, the 2.506 number is sqrt(2*pi)
      //.33f = the cutoff from sigma
      const float sigma_mul = 1.f; //.75f;
      const float r = distance / (sigma * sigma_mul);

      return (atomic_number / ( sigma * sigma_mul * 2.5066282746f )) * exp(-.5f * r * r);
    }
  };
  #endif
  
  struct ForcesRDF : public RadialDistanceFunction
  {
    virtual float evaluate(float distance_squared)
    {
      const float r2 = distance_squared * inv_sigma_squared;
      //cerr << "atomid= " << atomid << ", force mag = " << cd->forces[atomid].length3() << endl;
      return (cd->forces[atomid].length3() * expf(-r2));
    }
  };
  
 
  
  CubeData()
  {
    volume_dims = vec3i(0,0,0);
    ws2gs_scale = vec4(0,0,0,0);
    ws2gs_translate = vec4(0,0,0,0);
  
    voxels_per_macrocell = 0;
    voxels_per_angstrom = 4.f;
    bond_adjacencies = 0;
    volumes.reserve(4);
    rdf_pad_percent = 5;
    cvi = 0;
		
		write_gridded_particles = false;
	
	xyz_subtype = 0;
    
    bounds_min = vec3f(-9e999999f, -9e999999f, -9e999999f);
    bounds_max = vec3f(9e999999f, 9e999999f, 9e999999f);
    
    ws2gs_transform.identity();
    gs2ws_transform.identity();
    
    cs2ws_transform.identity();
    ws2cs_transform.identity();
    
    cs2gs_transform.identity();
    gs2cs_transform.identity();
  }
  ~CubeData(){}
  
  Array3Base* get_current_volume()
  {
    if (volumes.size() == 0)
      return 0;
    //cerr << "cvi = " << cvi << endl;
    //cerr << "get_current_volume() = " << size_t(volumes[cvi]) << endl;
    return volumes[cvi];
  }

  void parse_args(string str)
  {
	  int argc=1;
	  char cstr[4096];

	  strcpy(cstr, str.c_str());
	  cstr[strlen(cstr)-1] = 0;		//remove the trailing space.
	  for(int i=0; i<str.length(); i++)
		  if (cstr[i] == ' ') argc++;

	  //cerr << "argc = " << argc << endl;

	  char** argv = new char*[argc];
	  argv[0] = cstr;
	  int c=1;
	  for(int i=0; i<str.length(); i++)
	  {
		  if (cstr[i] == ' ')
		  { 
			  cstr[i] = 0;
			  argv[c] = cstr+i+1;
			  //cerr << "argv[" << c << "]=" << argv[c] << endl;
			  c++;
		  }
	  }

	parse_args(argc, argv);

	delete[] argv;
}
  
  void parse_args(int argc, char* argv[])
  {
    for(int i=1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        if (!strcmp(argv[i], "-atoms_string")) 
        {
          atoms_string = string(argv[++i]);
        }
        else if (!strcmp(argv[i], "-voxels_per_macrocell") || !strcmp(argv[i], "-vpm") || !strcmp(argv[i], "-mc_width") || !strcmp(argv[i], "-mcw")) 
        {
          voxels_per_macrocell = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-voxels_per_angstrom") || !strcmp(argv[i], "-vpa")) 
        {
          voxels_per_angstrom = atof(argv[++i]);
        }        
        else if (!strcmp(argv[i], "-ws2gs_scale") || !strcmp(argv[i], "-scale"))
        {
          ws2gs_scale[0] = atof(argv[++i]);
          ws2gs_scale[1] = atof(argv[++i]);
          ws2gs_scale[2] = atof(argv[++i]);
          ws2gs_scale[3] = 1.f;
        }
        else if (!strcmp(argv[i], "-bas_scale"))
        {
          ws2gs_scale[0] = 1.f/atof(argv[++i]);
          ws2gs_scale[1] = 1.f/atof(argv[++i]);
          ws2gs_scale[2] = 1.f/atof(argv[++i]);
          ws2gs_scale[3] = 1.f;
        }
        else if (!strcmp(argv[i], "-ws2gs_translate") || !strcmp(argv[i], "-translate"))
        {
          ws2gs_translate[0] = atof(argv[++i]);
          ws2gs_translate[1] = atof(argv[++i]);
          ws2gs_translate[2] = atof(argv[++i]);
          ws2gs_translate[3] = 0.f;          
        }
        else if (!strcmp(argv[i], "-bas_translate"))
        {
          ws2gs_translate[0] = -atof(argv[++i]) * ws2gs_scale[0];
          ws2gs_translate[1] = -atof(argv[++i]) * ws2gs_scale[1];
          ws2gs_translate[2] = -atof(argv[++i]) * ws2gs_scale[2];
          ws2gs_translate[3] = 0.f;          
        }
        else if (!strcmp(argv[i], "-rdf_pad_percent")) 
        {
          rdf_pad_percent = atof(argv[++i]);
        } 
        else if (!strcmp(argv[i], "-volume_dims"))
        {
          volume_dims[0] = atoi(argv[++i]);
          volume_dims[1] = atoi(argv[++i]);
          volume_dims[2] = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-read_forces"))
        {
          read_forces(string(argv[++i]));
        }
        else if (!strcmp(argv[i], "-xyz_subtype"))
        {
          xyz_subtype = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-bounds_min"))
        {
          bounds_min[0] = atof(argv[++i]);
          bounds_min[1] = atof(argv[++i]);
          bounds_min[2] = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-bounds_max"))
        {
          bounds_max[0] = atof(argv[++i]);
          bounds_max[1] = atof(argv[++i]);
          bounds_max[2] = atof(argv[++i]);
        }
      }
    }
  }
  
  void parse_args_after_read(int argc, char* argv[])
  {
    for(int i=1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        if (!strcmp(argv[i], "-create_rdf_volume"))
        {
          create_rdf_volume<RadialDistanceFunction>();
        }
        if (!strcmp(argv[i], "-create_forces_rdfv"))
        {
          //create_rdf_volume<RadialDistanceFunction>();
          create_rdf_volume<ForcesRDF>();
        }
        else if (!strcmp(argv[i], "-write_xyz"))
        {
          write_xyz();
        }
        else if (!strcmp(argv[i], "-write_volume"))
        {
          write_volume();
        }
        else if (!strcmp(argv[i], "-write_gridded_particles"))
        {
          write_gridded_particles = true;
        }
        else if (!strcmp(argv[i], "-set_current_volume") || !strcmp(argv[i], "-set_cvi"))
        {
          cvi = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-read_nhdr"))
        {
          if (i+1 < argc)
          {
            if (argv[i+1][0] != '-')
              read_nhdr(true, string(argv[++i]));
            else
              read_nhdr(true, string(filebase + "_rdf.nhdr"));
          }
          else
            read_nhdr(true, string(filebase + "_rdf.nhdr"));
        }
        else if (!strcmp(argv[i], "-create_uchar_volume"))
        {
          create_uchar_volume();
        }
        else if (!strcmp(argv[i], "-set_filebase"))
        {
          filebase = string(argv[++i]);
        }
      }
    }
  }  
  
  bool read_filename(string _filename)
  {
    filename = _filename;
    string fbstr = filename;
    //cerr_dbg(  "filename = " << filename << endl  ); 
    int pos = fbstr.find_last_of(".");
    string suffix;
    string prefix;
    if (pos > 0)
    {
      suffix = fbstr.substr(pos);
      if (suffix == ".raw" || suffix == ".xyz" || suffix == ".nhdr" || suffix == ".hdr" || suffix == ".cube" || suffix == ".chgcar" || suffix == ".dat" )
      {
        prefix = fbstr.substr(0, pos);
        filebase = prefix.c_str();
        cerr_dbg(  "Stripping extension " << suffix << endl  ); 
      }
      else
      {
        cerr << "Unrecognized file format " << suffix << endl;
        return false;
      }
      
    }
    else
    {
      filebase = filename;
    }
    
    cerr_dbg(  "filebase = " << filebase << endl  ); 
    
    bool result;
    
    if (suffix == ".nhdr")
    {
      result = read_nhdr();
    }
    else if (suffix == ".xyz")
    {
      result = read_xyz();
    }    
    else if (suffix == ".cube")
    {
      result = read_cube();
    }
    else if (suffix == ".chgcar" || filename == "CHGCAR")
    {
      result = read_chgcar();
    }
    else if (suffix == ".dat")
    {
      result = read_dat();
    }
    else
    {
      cerr_dbg(  "CubeData not handled!" << endl  ); 
      return false;
    }
    
    if (result == false)
    {
      cerr << "Miscellaneous error reading data from " << filename << endl;
      return false;
    }
    
    return true;
  }
 
  
  bool read_nhdr(bool post_read_nhdr = false, string nhdr_filename = "")
  {
    cerr << "read_nhdr" << endl;

    if (nhdr_filename == "")
      nhdr_filename = string(filebase + ".nhdr");

    FILE* fin = fopen(nhdr_filename.c_str(), "r+");
    
    if (!fin)
    {
      cerr << "Failed to load " << nhdr_filename << endl;
      exit(0);
    }
    
    char buffer[1024];
    char buffer2[1024];
    char buffer3[1024];
    char buffer4[1024];
    float x;
    
    bool nrrd_uchar_data = false;

    string raw_filename;
    string xyz_filename;
    
    if (post_read_nhdr)
    {
      ws2gs_scale = vec4(0,0,0,0);
      ws2gs_translate = vec4(0,0,0,0);
    }

    Array3Base a3b;    
    a3b.range = vec2f(0,0);
    
    while(fgets(buffer, 1024, fin))
    {
      if (strstr(buffer, "sizes:"))
      {
        sscanf(buffer, "%s %d %d %d", buffer2, &volume_dims[0], &volume_dims[1], &volume_dims[2]); 
      }
      else if (strstr(buffer, "#-"))
      {
        char* bp = &buffer[0]+1;
        buffer[strlen(buffer)-1] = ' ';	//replace newline with a space
        cerr << "nhdr options string: " << bp << endl;
        nhdr_options_string += string(bp);
      }
      else if (strstr(buffer, "type:"))
      {
        if (strstr(buffer, "unsigned char") || strstr(buffer, "uchar"))
        {
          a3b.type = string("uchar");
          nrrd_uchar_data = true;
        }
        else
        {
          a3b.type = string("float");
          nrrd_uchar_data = false;          
        }
      }
      else if (strstr(buffer, "min:"))
      {
        sscanf(buffer, "%s %f", buffer2, &x);
        a3b.range[0] = x;
      }
      else if (strstr(buffer, "max:"))
      {
        sscanf(buffer, "%s %f", buffer3, &x);
        a3b.range[1] = x;
      }
      else if (strstr(buffer, "old min:"))
      {
        sscanf(buffer, "%s %s %f", buffer2, buffer3, &x);
        a3b.old_range[0] = x;
      }
      else if (strstr(buffer, "old max:"))
      {
        sscanf(buffer, "%s %s %f", buffer2, buffer3, &x);
        a3b.old_range[1] = x;
      }         
      else if (strstr(buffer, "data file:"))
      {
        sscanf(buffer, "%s %s %s", buffer2, buffer3, buffer4);
        //string filename = string(buffer4);
        //a3b.filebase = filename.substr(0, filename.find_last_of("."));
	a3b.filebase = filebase.substr(0, filename.find_last_of("."));
      }   
      else if (strstr(buffer, "#xyz file:"))
      {
        sscanf(buffer, "%s %s %s", buffer2, buffer3, buffer4);
        xyz_filename = string(buffer4);
      }
      else if (strstr(buffer, "#ws2gs_translate:"))
      {
        sscanf(buffer, "%s %f %f %f\n", buffer2, &ws2gs_translate[0], &ws2gs_translate[1], &ws2gs_translate[2]);
      }
      else if (strstr(buffer, "#ws2gs_scale:"))
      {
        sscanf(buffer, "%s %f %f %f\n", buffer2, &ws2gs_scale[0], &ws2gs_scale[1], &ws2gs_scale[2]);
      }
    }
    fclose(fin);
    
    if (post_read_nhdr)
    {
      if (ws2gs_scale[1] == 0.f)
         ws2gs_scale[1] = ws2gs_scale[0];
      if (ws2gs_scale[2] == 0.f)
	 ws2gs_scale[2] = ws2gs_scale[0];
      ws2gs_scale[3] = 1.f;
      ws2gs_translate[3] = 0.f;
    }
                
    cerr << "volume_dims = " << volume_dims << endl;
    cerr << "volume range = " << a3b.range << endl;
    cerr << "volume old range = " << a3b.old_range << endl;
    cerr << "ws2gs_translate = " << ws2gs_translate << endl;
    cerr << "ws2gs_scale = " << ws2gs_scale << endl;
    cerr << "xyz_filename = " << xyz_filename << endl;

    parse_args(nhdr_options_string);

    nhdr_options_string = "";
    
    a3b.dims = volume_dims;
    Array3Base* volume;
    if (nrrd_uchar_data)
      volume = new Array3<unsigned char>(&a3b);
    else  
      volume = new Array3<float>(&a3b);
      
    volume->read_raw();
    cus->display_uchar_volume = nrrd_uchar_data;

    volumes.push_back(volume); 
   
    if (volume->range[0] >= volume->range[1])
    {
      cerr << "finding min/max" << endl;
      volume->compute_minmax();
      cerr << "volume range = " << volume->range << endl;
    }

    if (post_read_nhdr || xyz_filename != "")
    {
      volume->xyz_filebase = xyz_filename.substr(0, xyz_filename.find_last_of("."));
      read_xyz(xyz_filename);
    }
    
    return true;
  }

  bool read_eff(string _filename = "")
  {
    if (_filename == "")
      _filename = filebase + ".eff";
    
    FILE* file = fopen(_filename.c_str(), "r");
    
    if (!file)
    {
      cerr << "Couldn't find " << filename << endl;
      return false;
    }
    
    int npts;

    
    fscanf(file,"%d",&npts);

    char idbuf[1024];
    
    if (xyz_subtype)
      fgets(idbuf, 1024, file);
      
    cerr << "Reading xyz molecule with " << npts << " vertices." << endl;

    if (cus)
    {
      if (cus->max_atoms && cus->max_atoms < npts)
      {
        cerr << "But, user told us to use only the first " << cus->max_atoms << " atoms." << endl; 
        npts = cus->max_atoms;
      }
    }
    
    ws_atoms.resize(npts);
 
 
    char identifier[512];
    float x,y,z,w;            
    
    fgets(idbuf, 1024, file);
    cerr << "idbuf = " << idbuf << endl;

		sort(special_atoms.begin(), special_atoms.end());
		cerr << "special_atoms.size() = " << special_atoms.size() << endl;
    
		int sai=0;
    int points_read = 0;
    
    vec4f min_atom = vec4f(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    vec4f max_atom = vec4f(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    
    float tmpE0, tmpE1;
    int tmpID;

    for(int i = 0; i < ws_atoms.size(); i++)
    {
      if (xyz_subtype)
      {
        fscanf(file,"%s %f %f %f %d %f %f\n", identifier, &x,&y,&z,&tmpID,&tmpE0,&tmpE1);
      }

      char tmpid[4];
      if (cus->use_custom_atoms)
      {
        sprintf(tmpid, "_%s", identifier);
        strcpy(identifier, tmpid);
      }
      else
        fscanf(file,"%s %f %f %f\n", identifier, &x,&y,&z);
			if (sai < special_atoms.size() && i == special_atoms[sai])
			{
				sai++;
				strcpy(identifier, "SP");
			}
      
      w = cus->get_atomic_number(string(identifier));
      
      if (w >= cus->abort_atomic_number)
      {
        return false;
      }
      
      //cout << x << " " << y << " " << z << " " << w << endl;

      if (x > bounds_min[0] && x < bounds_max[0] && y > bounds_min[1] && y < bounds_max[1] && z > bounds_min[2] && z < bounds_max[2])
	  {
        ws_atoms[i] = vec4f(x,y,z,w);
        min_atom = min(ws_atoms[i], min_atom);
        max_atom = max(ws_atoms[i], max_atom);
	  }
      else
        cerr << "out of range atom: " << vec4f(x,y,z,w) << endl;
      //Matrix44::transform(gs_atoms[i], ws2gs_transform, ws_atoms[i]);
      
      points_read++;
    }
          
    cerr << "points read = " << points_read << " / " << ws_atoms.size() << endl;
    cerr << "done!" << endl;    
      
    fclose(file);

    //NOTE, even if we don't have volume data, we still want a "grid" scale from which to make macrocells
    //this is determined by command-line options

	  //if ws2gs_scale != 0, then user specifies a transform and dims. Else, we need to create one
    if (ws2gs_scale[0] == 0.f)
    {
      cerr << "Warning, read_xyz() was called with no volume or transform specified. Generating one." << endl;
      
      vec4 pad_vector = (max_atom - min_atom) * (rdf_pad_percent * .01f * .5f);
      cerr << "min_atom = " << min_atom << endl;
      cerr << "max_atom = " << max_atom << endl;
      cerr << "rdf_pad_percent = " << rdf_pad_percent << endl;
      vec4 max_atom_pad = max_atom + pad_vector;
      vec4 min_atom_pad = min_atom - pad_vector;
      
      vec4 pad_extents = max_atom_pad - min_atom_pad;
      cerr << "pad_extents = " << pad_extents << endl;
      
      if (volume_dims[0])
      {
        voxels_per_angstrom = volume_dims[0] / pad_extents[0];
        ws2gs_scale[0] = voxels_per_angstrom;
      }

      cerr << "voxels_per_angstrom = " << voxels_per_angstrom << endl;
      //if volume_dims == 0, we need to determine them
      volume_dims[0] = int(ceilf(pad_extents[0] * voxels_per_angstrom));
      volume_dims[1] = int(ceilf(pad_extents[1] * voxels_per_angstrom));
      volume_dims[2] = int(ceilf(pad_extents[2] * voxels_per_angstrom));
      cerr << "volume_dims = " << volume_dims << endl;
      
      ws2gs_scale = vec4(voxels_per_angstrom);
      ws2gs_scale[3] = 1.f;
      
      ws2gs_translate = -min_atom_pad * ws2gs_scale;
      ws2gs_translate[3] = 0.f;
      
    }
    else //ws2gs_scale and ws2gs_translate have been specified
    {
    
      if (ws2gs_scale[1] == 0.f)
			  ws2gs_scale[1] = ws2gs_scale[0];
		  if (ws2gs_scale[2] == 0.f)
			  ws2gs_scale[2] = ws2gs_scale[0];
			  
			voxels_per_angstrom = ws2gs_scale[0];
    }
			
    cerr << "ws2gs_translate = " << ws2gs_translate << endl;
    cerr << "ws2gs_scale = " << ws2gs_scale << endl;
    
    ws2gs_transform.identity();
    ws2gs_transform.m[0][0] = ws2gs_scale[0];
    ws2gs_transform.m[1][1] = ws2gs_scale[1];
    ws2gs_transform.m[2][2] = ws2gs_scale[2];
    ws2gs_transform.m[0][3] = ws2gs_translate[0];
    ws2gs_transform.m[1][3] = ws2gs_translate[1];
    ws2gs_transform.m[2][3] = ws2gs_translate[2];
    ws2gs_transform.scrub();
    
    cerr << "ws2gs_transform = " << ws2gs_transform << endl;
    
    Matrix44::inverse(gs2ws_transform, ws2gs_transform);

    cerr << "gs2ws_transform = " << gs2ws_transform << endl;
   
    cerr << "Creating gs_atoms..." << endl;
    int sqrt_size = int(ceilf(sqrt(ws_atoms.size())));
#ifdef BAS_TEXTURE_BUFFER
    gs_atoms.reserve(ws_atoms.size());
#else
    gs_atoms.reserve(sqrt_size * sqrt_size);
#endif
    gs_atoms.resize(ws_atoms.size());
          
    for(int i=0; i<ws_atoms.size(); i++)
    {
      gs_atoms[i] = transform_atom(ws2gs_transform, ws_atoms[i]);
      //gs_atoms[i] += vec4(.5f, .5f, .5f, 0.f);
      
      //if (i < 20)
      //  cerr_dbg( "gs_atoms[" << i << "] = " << gs_atoms[i] << endl );
    }
    
    cerr << "done." << endl;
    
    return true;
  }
    
  bool read_xyz(string _filename = "")  
  {
    if (_filename == "")
      _filename = filebase + ".xyz";
	  
	cerr << "xyz_subtype = " << xyz_subtype << endl;
    
    FILE* file = fopen(_filename.c_str(), "r");
    
    if (!file)
    {
      cerr << "Couldn't find " << filename << endl;
      return false;
    }
    
    int npts;
    char idbuf[1024];

    if (xyz_subtype==2)
    {
      fgets(idbuf, 1024, file);
      fgets(idbuf, 1024, file);
      fgets(idbuf, 1024, file);
    }
  
    if (xyz_subtype == 3)
    {
      fscanf(file,"%d %s",&npts, idbuf);
      fgets(idbuf, 1024, file);
    }  
    else
      fscanf(file,"%d",&npts);

    
    if (xyz_subtype<=1)
      fgets(idbuf, 1024, file);
    else if (xyz_subtype==2)
    {
      for(int i=0; i<5; i++)
        fgets(idbuf, 1024, file);   
    }
      
    cerr << "Reading xyz molecule with " << npts << " vertices." << endl;

    if (cus)
    {
      if (cus->max_atoms && cus->max_atoms < npts)
      {
        cerr << "But, user told us to use only the first " << cus->max_atoms << " atoms." << endl; 
        npts = cus->max_atoms;
      }
    }
    
    ws_atoms.resize(npts);
 
 
    char identifier[512];
    float x,y,z,w;            

    fgets(idbuf, 1024, file);
    cerr << "idbuf = " << idbuf << endl;
    
#if 0
    //AARONBAD -- this should NOT be in read_xyz. This should be in the nhdr as a comment
    char* transs = strstr(idbuf, "translate");
    if (transs)
      sscanf(transs, "%s %f %f %f", identifier, &ws2gs_translate[0], &ws2gs_translate[1], &ws2gs_translate[2]);
    
    char* scales = strstr(idbuf, "scale");
    if (scales)
      sscanf(scales, "%s %f %f %f", identifier, &ws2gs_scale[0], &ws2gs_scale[1], &ws2gs_scale[2]);
#endif

  
		sort(special_atoms.begin(), special_atoms.end());
		cerr << "special_atoms.size() = " << special_atoms.size() << endl;
    
		int sai=0;
    int points_read = 0;
    
    vec4 min_atom = vec4(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    vec4 max_atom = vec4(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
    
    float tmpE0, tmpE1;
    int tmpID;

    for(int i = 0; i < ws_atoms.size(); i++)
    {
      if (xyz_subtype==1 || xyz_subtype==3)
      {
        fscanf(file,"%s %f %f %f %d %f %f\n", identifier, &x,&y,&z,&tmpID,&tmpE0,&tmpE1);
      }

      char tmpid[4];
      if (cus->use_custom_atoms)
      {
        sprintf(tmpid, "_%s", identifier);
        strcpy(identifier, tmpid);
      }
      else
        fscanf(file,"%s %f %f %f\n", identifier, &x,&y,&z);
      if (xyz_subtype==2)
        strcpy(identifier, "C");
			if (sai < special_atoms.size() && i == special_atoms[sai])
			{
				sai++;
				strcpy(identifier, "SP");
			}
      
      w = cus->get_atomic_number(string(identifier));
      
      if (w >= cus->abort_atomic_number)
      {
        return false;
      }
      
      //cout << x << " " << y << " " << z << " " << w << endl;

      if (x > bounds_min[0] && x < bounds_max[0] && y > bounds_min[1] && y < bounds_max[1] && z > bounds_min[2] && z < bounds_max[2])
	  {
        ws_atoms[i] = vec4f(x,y,z,w);
      	min_atom = min(ws_atoms[i], min_atom);
      	max_atom = max(ws_atoms[i], max_atom);
	  }
      else
        cerr << "out of range atom: " << vec4f(x,y,z,w) << endl;
      //Matrix44::transform(gs_atoms[i], ws2gs_transform, ws_atoms[i]);
      

      points_read++;
    }
          
    cerr << "points read = " << points_read << " / " << ws_atoms.size() << endl;
    cerr << "done!" << endl;    
      
    fclose(file);

    //NOTE, even if we don't have volume data, we still want a "grid" scale from which to make macrocells
    //this is determined by command-line options

	  //if ws2gs_scale != 0, then user specifies a transform and dims. Else, we need to create one
    if (ws2gs_scale[0] == 0.f)
    {
      cerr << "Warning, read_xyz() was called with no volume or transform specified. Generating one." << endl;
      
      vec4 pad_vector = (max_atom - min_atom) * (rdf_pad_percent * .01f * .5f);
      cerr << "min_atom = " << min_atom << endl;
      cerr << "max_atom = " << max_atom << endl;
      cerr << "rdf_pad_percent = " << rdf_pad_percent << endl;
      vec4 max_atom_pad = max_atom + pad_vector;
      vec4 min_atom_pad = min_atom - pad_vector;
      
      vec4 pad_extents = max_atom_pad - min_atom_pad;
      cerr << "pad_extents = " << pad_extents << endl;
      
      if (volume_dims[0])
      {
        voxels_per_angstrom = volume_dims[0] / pad_extents[0];
        ws2gs_scale[0] = voxels_per_angstrom;
      }

      cerr << "voxels_per_angstrom = " << voxels_per_angstrom << endl;
      //if volume_dims == 0, we need to determine them
      volume_dims[0] = int(ceilf(pad_extents[0] * voxels_per_angstrom));
      volume_dims[1] = int(ceilf(pad_extents[1] * voxels_per_angstrom));
      volume_dims[2] = int(ceilf(pad_extents[2] * voxels_per_angstrom));
      cerr << "volume_dims = " << volume_dims << endl;
      
      ws2gs_scale = vec4(voxels_per_angstrom);
      ws2gs_scale[3] = 1.f;
      
      ws2gs_translate = -min_atom_pad * ws2gs_scale;
      ws2gs_translate[3] = 0.f;
      
    }
    else //ws2gs_scale and ws2gs_translate have been specified
    {
    
      if (ws2gs_scale[1] == 0.f)
			  ws2gs_scale[1] = ws2gs_scale[0];
		  if (ws2gs_scale[2] == 0.f)
			  ws2gs_scale[2] = ws2gs_scale[0];
			  
			voxels_per_angstrom = ws2gs_scale[0];
    }
			
    cerr << "ws2gs_translate = " << ws2gs_translate << endl;
    cerr << "ws2gs_scale = " << ws2gs_scale << endl;
    
    ws2gs_transform.identity();
    ws2gs_transform.m[0][0] = ws2gs_scale[0];
    ws2gs_transform.m[1][1] = ws2gs_scale[1];
    ws2gs_transform.m[2][2] = ws2gs_scale[2];
    ws2gs_transform.m[0][3] = ws2gs_translate[0];
    ws2gs_transform.m[1][3] = ws2gs_translate[1];
    ws2gs_transform.m[2][3] = ws2gs_translate[2];
    ws2gs_transform.scrub();
    
    cerr << "ws2gs_transform = " << ws2gs_transform << endl;
    
    Matrix44::inverse(gs2ws_transform, ws2gs_transform);

    cerr << "gs2ws_transform = " << gs2ws_transform << endl;
   
    cerr << "Creating gs_atoms..." << endl;
    int sqrt_size = int(ceilf(sqrt(ws_atoms.size())));
#ifdef BAS_TEXTURE_BUFFER
    gs_atoms.reserve(ws_atoms.size());
#else
    gs_atoms.reserve(sqrt_size * sqrt_size);
#endif
    gs_atoms.resize(ws_atoms.size());
          
    for(int i=0; i<ws_atoms.size(); i++)
    {
      gs_atoms[i] = transform_atom(ws2gs_transform, ws_atoms[i]);
      //gs_atoms[i] += vec4(.5f, .5f, .5f, 0.f);
      
      //if (i < 20)
      //  cerr_dbg( "gs_atoms[" << i << "] = " << gs_atoms[i] << endl );
    }
    
    cerr << "done." << endl;
    
    return true;
  }
    
  bool read_chgcar(string filename = "")
  {
    if (filename == "")
      filename = filebase + ".chgcar";
    else
      filebase = filename.substr(0, filename.find_last_of(".")); 
    
    FILE* fin = fopen(filename.c_str(), "r");
    
    if (!fin)
    {
      printf("Couldn't open %s, trying with .chgcar", filename.c_str());
      
      filename = filebase + ".chgcar";
      fin = fopen(filename.c_str(), "r");
      
      if (!fin)
      {
        printf("Nope! Exiting...");
        return 0;   
      }
      
      printf("read %s", filename.c_str());
    }
     
    const int BSIZE = 1024;
    char buffer[BSIZE];
    fgets(buffer, BSIZE, fin);
    for(char* c= buffer; *c != 0; c++)
    {
      if (*c == '\n') 
        *c = ' ';
    }
    description = string(buffer);
    
    float charge_multiplier;
    fscanf(fin, "%f \n", &charge_multiplier);
    
    cs2ws_transform.identity();
    fscanf(fin, "%f %f %f \n", &cs2ws_transform.m[0][0], &cs2ws_transform.m[1][0], &cs2ws_transform.m[2][0]);
    fscanf(fin, "%f %f %f \n", &cs2ws_transform.m[0][1], &cs2ws_transform.m[1][1], &cs2ws_transform.m[2][1]);
    fscanf(fin, "%f %f %f \n", &cs2ws_transform.m[0][2], &cs2ws_transform.m[1][2], &cs2ws_transform.m[2][2]);
    cs2ws_transform.m[3][3] = 1.f;
    
    cerr << "cs2ws_transform = " << cs2ws_transform << endl;
    
    int num_atoms = 0;
    
    //fscanf(fin, "%d", &num_atoms);

    //fscanf(fin, "%s", buffer);
    char buffer2[BSIZE];
    fgets(buffer, BSIZE, fin);
    sscanf(buffer, "%s", buffer2);
    //cerr_dbg(  "buffer = " << buffer << endl  ); 
    if (atoi(buffer2) > 0) //if we just have numbers, no atoms
    {
      cerr_dbg(  "atoi(buffer2) = " << atoi(buffer2) << endl  ); 
      //cerr << "atoms_string = " << atoms_string << endl;
      if (atoms_string == "")
      {
        cerr <<  "Fatal: cubeutils requires atom types in CHGCAR. Put these below the unit cell matrix, or use the argument -atoms_string [atom0] [atom1] [atom2] etc." << endl; 
        return false;
      }
    }
    else
    {
      atoms_string = string(buffer);
      fgets(buffer, BSIZE, fin);
    }
    
    cerr_dbg(  "atoms_string = " << atoms_string << endl  ); 
    
    //parse the atom type string
    {
      stringstream sstr(atoms_string);
      string word;
      while(getline(sstr, word, ' '))
      {
        if (strcmp(word.c_str(), ""))
        {
          char buf[16];
          sscanf(word.c_str(), "%s", buf);
          cerr_dbg(  "atom type = " << buf << endl  ); 
          cerr_dbg(  "atomic# = " << cus->get_atomic_number(buf) << endl  ); 
          atom_type_counts.push_back(vec2i(cus->get_atomic_number(buf), 0));
        }
      }
    }
    
    //now read the number of each atom
    {
      stringstream sstr(buffer);
      string word;
      int i=0;
      while(getline(sstr, word, ' '))
      {
        if (strcmp(word.c_str(), ""))
        {
          cerr_dbg(  "atom count = " << word << endl  ); 
          //atom_type_counts[i][0] = 0;
          atom_type_counts[i][1] = atoi(word.c_str());
          num_atoms += atom_type_counts[i][1];
          i++;
        }
      }
    }
    
    cerr_dbg(  "num_atoms = " << num_atoms << endl  ); 
        
    fscanf(fin, "%s", buffer);  //should say something like "Direct", just read it in and ignore

    cs_atoms.reserve(num_atoms);
    ws_atoms.reserve(num_atoms);
    
    for(int c=0; c<atom_type_counts.size(); c++)
    {
      for(int i=0; i<atom_type_counts[c][1]; i++)
      {
        //in chgcar, atoms are in unit cell space by default
        vec4f atom;
        fscanf(fin, "%f %f %f\n", &atom[0], &atom[1], &atom[2]);
        
        atom[3] = float(atom_type_counts[c][0]);
        //cerr_dbg(  "atom[3] = " << atom[3] << endl  );
        cs_atoms.push_back(atom);
        vec4 wsatom = atom;
        wsatom[3] = 0.f;
        wsatom = cs2ws_transform * wsatom;
        wsatom[3] = atom[3];
        ws_atoms.push_back(wsatom);
        //cerr <<  "atoms[" << atoms.size()-1 << "] = " << atoms[atoms.size()-1] << endl;
      }
    }
        
    //skip a line
    //fgets(buffer, BSIZE, fin);
    
    fscanf(fin, "%d %d %d\n", &volume_dims[0], &volume_dims[1], &volume_dims[2]);
    cerr_dbg(  "volume_dims = " << volume_dims << endl  ); 
    
    
    Matrix44::inverse(ws2cs_transform, cs2ws_transform);
    
    //so we have cs2ws and ws2cs
    //we need ws2gs and cs2gs
    //NOTE: we use cs2gs to build gs_atoms
    
    cs2gs_transform.identity();
    cs2gs_transform.scale( vec4(volume_dims[0], volume_dims[1], volume_dims[2], 1.f) );
    Matrix44::inverse(gs2cs_transform, cs2gs_transform);
    cerr_dbg(  "cs2gs_transform = " << cs2gs_transform << endl  ); 
    
    ws2gs_transform = ws2cs_transform;
    ws2gs_transform.scale( vec4(volume_dims[0], volume_dims[1], volume_dims[2], 1.f) );
    Matrix44::inverse(gs2ws_transform, ws2gs_transform);
    cerr_dbg(  "ws2gs_transform = " << ws2gs_transform << endl  ); 
    
    const int num_voxels = volume_dims[0] * volume_dims[1] * volume_dims[2];
    
    Array3<float>* volume = new Array3<float>;
    volumes.push_back(volume);
    volume->dims = volume_dims;
    volume->type = "float";
    volume->filebase = filebase;
    
    volume->resize(volume_dims[0], volume_dims[1], volume_dims[2]);
        
    //divide value by cell volume to get rho -- see http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.1151
    //AARONBAD -- do we want the sqrt here? Maria said yes, but web says no.
    
    double cell_volume;
    
    
    //cell_volume = double(ws2gs_transform.m[0][0]) * double(ws2gs_transform.m[1][1]) * double(ws2gs_transform.m[2][2]) / (double(cs2gs_transform.m[0][0]) * double(cs2gs_transform.m[1][1]) * double(cs2gs_transform.m[2][2]));
    //cell_volume = double(volume_dims[0] * volume_dims[1] * volume_dims[2]) / (double(cs2ws_transform.m[0][0]) * double(cs2ws_transform.m[1][1]) * double(cs2ws_transform.m[2][2]));
    //cell_volume /= double(charge_multiplier);
    
    //if (cus->use_sqrt)
    //  cell_volume = sqrt(cell_volume);

    //charges are per voxel; we need to convert to per Ang^3
    float scale_cell_volume = charge_multiplier * ws2cs_transform.m[0][0] * ws2cs_transform.m[1][1] * ws2cs_transform.m[2][2];
    
    if (cus->use_sqrt)
      scale_cell_volume = sqrtf(scale_cell_volume);
    
    volume->range[0] = std::numeric_limits<float>::max();
    volume->range[1] = std::numeric_limits<float>::min();
    
    for(int i=0; i<num_voxels; i++)
    {
      float val;
      fscanf(fin, "%f", &val);
      
      //printf("%f ", val); 
      
      val *= scale_cell_volume; 
      
      //val = max(val, 0.f);
      if (cus->use_sqrt)
      {
        val = max(val, 0.f);
        val = sqrtf(val); //convert from charge density (rho) to wavefunction (sqrt(rho))
      }

      volume->range[0] = min(volume->range[0], val);
      volume->range[1] = max(volume->range[1], val);
      volume->data[i] = val;
    }
    
    fclose(fin);
        
    printf("volume->range is [%f, %f]\n", volume->range[0], volume->range[1]);
    

    
    
    //note, we want to build gs_atoms manually, because we do NOT want the half-cell offset
    gs_atoms.resize(ws_atoms.size());
    for(int i=0; i<ws_atoms.size(); i++)
    {
      gs_atoms[i] = transform_atom(ws2gs_transform, ws_atoms[i]);
      //if (i < 20)
      //  cerr_dbg( "gs_atoms[" << i << "] = " << gs_atoms[i] << endl );
    }
        
    return true;
  }

  bool read_cube(string filename = "")
  {
    if (filename == "")
      filename = filebase + ".cube";
    else
      filebase = filename.substr(0, filename.find_last_of("."));
      

    FILE* fin = fopen(filename.c_str(), "r");
    if (!fin)
    {
      printf("Couldn't open %s", filename.c_str());
      return 0;   
    }
    
    const int BSIZE = 512;
    char buffer[BSIZE];

    fgets(buffer, BSIZE, fin);
    for(char* c= buffer; *c != 0; c++)
    {
      if (*c == '\n') 
        *c = ' ';
    }    
    name = string(buffer);
    cerr_dbg(  "name = " << name << endl  ); 

    fgets(buffer, BSIZE, fin);
    for(char* c= buffer; *c != 0; c++)
    {
      if (*c == '\n') 
        *c = ' ';
    }    
    description = string(buffer);
    cerr_dbg(  "description = " << description << endl  ); 
    
    int num_atoms;
    vec4 translate;
    fscanf(fin, "%d %f %f %f\n", &num_atoms, &translate[0], &translate[1], &translate[2]);
    fscanf(fin, "%d %f %f %f\n", &volume_dims[0], &gs2ws_transform.m[0][0], &gs2ws_transform.m[0][1], &gs2ws_transform.m[0][2]);
    fscanf(fin, "%d %f %f %f\n", &volume_dims[1], &gs2ws_transform.m[1][0], &gs2ws_transform.m[1][1], &gs2ws_transform.m[1][2]);
    fscanf(fin, "%d %f %f %f\n", &volume_dims[2], &gs2ws_transform.m[2][0], &gs2ws_transform.m[2][1], &gs2ws_transform.m[2][2]);
    
    gs2ws_transform.m[0][3] = gs2ws_transform.m[1][3] = gs2ws_transform.m[2][3] = 0.f;
    gs2ws_transform.m[0][3] = translate[0];
    gs2ws_transform.m[1][3] = translate[1];
    gs2ws_transform.m[2][3] = translate[2];
    gs2ws_transform.m[3][3] = 1.f;
    
    cerr_dbg(  "gs2ws_transform = " << gs2ws_transform << endl  ); 
    
    Matrix44::inverse(ws2gs_transform, gs2ws_transform);
    
    cerr_dbg(  "ws2gs_transform = " << ws2gs_transform << endl  ); 
    
    cerr_dbg(  "num_atoms = " << num_atoms << endl  ); 
    
    ws_atoms.reserve(num_atoms);
    for(int i=0; i<num_atoms; i++)
    {
      //in cube, atoms are in world space by default      
      vec4f atom;
      float tmp;
      fscanf(fin, "%f %f %f %f %f\n", &atom[3], &tmp, &atom[0], &atom[1], &atom[2]);
      cerr_dbg(  "atom = " << atom << endl  ); 
      ws_atoms.push_back(atom);
      //vec4 gs_atom = ws2gs_transform * vec4(atom[0], atom[1], atom[2], 1.f);
      //cerr_dbg(  "gs_atom = " << gs_atom << endl  ); 
      //gs_atoms.push_back(gs_atom);
    }
    
    Array3<float>* volume = new Array3<float>();
    volumes.push_back(volume);
    volume->dims = volume_dims;
    volume->type = "float";
    volume->filebase = filebase;
    
    const int num_voxels = volume_dims[0] * volume_dims[1] * volume_dims[2];
    volume->resize(volume_dims[0], volume_dims[1], volume_dims[2]);

    volume->range[0] = std::numeric_limits<float>::max();
    volume->range[1] = std::numeric_limits<float>::min();
    
    //Note: .cube's are written in z-fast order instead of x-fast.
    for(int i=0; i<volume_dims[0]; i++)
      for(int j=0; j<volume_dims[1]; j++)
        for(int k=0; k<volume_dims[2]; k++)
        {
          float val;
          fscanf(fin, "%f", &val);
          
          if (cus->use_sqrt)
          {
            val = max(val, 0.f);
            val = sqrtf(val); //convert from charge density (rho) to wavefunction (sqrt(rho))
          }
          
          volume->range[0] = min(volume->range[0], val);
          volume->range[1] = max(volume->range[1], val);
          volume->get_data(i,j,k) = val;
          
        }

    cerr << "volume range = " << volume->range << endl; 
    
    //note, we want to build gs_atoms manually, because we do NOT want the half-cell offset
    gs_atoms.resize(ws_atoms.size());
    for(int i=0; i<ws_atoms.size(); i++)
    {
      gs_atoms[i] = transform_atom(ws2gs_transform, ws_atoms[i]);
      //if (i < 20)
      //  cerr_dbg( "gs_atoms[" << i << "] = " << gs_atoms[i] << endl );
    }
        
    return true;
  }
  
  bool read_forces(string _filename)  
  {
    forces_filename = _filename;
    FILE* file = fopen(forces_filename.c_str(), "r");
    
    if (!file)
    {
      cerr << "Couldn't find " << forces_filename << endl;
      return false;
    }
    
    int npts;
    int points_read = 0;
    
    fscanf(file,"%d",&npts);
    cerr << "Reading forces file with " << npts << " vertices." << endl;
    
    char identifier[512];
    float x,y,z,w;            
    fgets(identifier, 512, file);
    
    forces.resize(npts);
    for(int i = 0; i < forces.size(); i++)
    {
			fscanf(file,"%s %f %f %f", identifier, &x,&y,&z);
      w = cus->get_atomic_number(string(identifier));
      forces[i] = vec4f(x,y,z,w);
      cerr << "forces[" << i << "] = " << forces[i] << endl;
      points_read++;
    }
          
    cerr << "points read = " << points_read << " / " << forces.size() << endl;
    cerr << "done!" << endl;    
      
    fclose(file);

    return true;
  }

#if 1
  template<typename RDF>
  void create_rdf_volume()
  {
    cerr << "create_rdf_volume" << endl;
    //creates an approximate charge density volume
    
    cerr << "ws2gs_transform = " << ws2gs_transform;
    cerr << "volume_dims = " << volume_dims << endl;
    
   
    Array3<float>* rdf_volume = new Array3<float>;
    rdf_volume->dims = volume_dims;
    rdf_volume->type = "float";
    rdf_volume->filebase = filebase + "_rdf";
    rdf_volume->allocate();
    rdf_volume->initialize(0);
    volumes.push_back(rdf_volume);
    
    const vec3 vdm1_vec3 = volume_dims - vec3(1,1,1);
    const vec4 vdm1 = vec4(vdm1_vec3[0], vdm1_vec3[1], vdm1_vec3[2], 16384);
    
    Clock clock;
    clock.start();
    
    //if (omp_get_dynamic())
    //omp_set_dynamic(1);
    
    //omp_set_dynamic(0);
    //int ncores=omp_get_max_threads();  

    //omp_set_num_threads( ncores );
    //cerr << "ncores = " << ncores << endl;
    //cerr << "omp_get_num_threads = " << omp_get_num_threads() << endl;
  
#pragma omp parallel for
    for(int h=0; h<gs_atoms.size(); h++)
    {
      vec4 p = ws_atoms[h];
      
      //if (h % 1000 == 0)
      //  cerr << "h = " << h << endl;
      
      RDF rdf;
      rdf.init(this, cus, h, int(p[3]));
      
      if (rdf.atom == 255)
        continue;
        
      //cerr << "ws_atoms[" << h << "] = " << p << endl;
      
      //find the min, max extents of this point's contribution
      const float clamp_radius = rdf.get_clamp_radius();
      //cerr << "clamp_radius = " << clamp_radius << endl;
      vec4 cr4 = vec4(clamp_radius, clamp_radius, clamp_radius, 0.f);
      
      //cerr << "clamp_radius = " << clamp_radius << endl;
      
      vec4 emin = p - cr4;
      vec4 emax = p + cr4;
      emin[3] = emax[3] = 1.f;
      
      //cerr << "emin = " << emin << endl;
      //cerr << "emax = " << emax << endl;
      
      //ok now these extents are in Angstroms. We need them in voxel space
      vec4 gs_emin = ws2gs_transform * emin;
      vec4 gs_emax = ws2gs_transform * emax;

      gs_emin = max(gs_emin, vec4(0,0,0,0));
      gs_emax = min(gs_emax, vdm1);
    
      vec4f v;
      for(int k=gs_emin[2]; k<gs_emax[2]; k++)
      {
        v[2] = float(k);
        for(int j=gs_emin[1]; j<gs_emax[1]; j++)
        {
          v[1] = float(j);
          for(int i=gs_emin[0]; i<gs_emax[0]; i++)
          {
             v[0] = float(i);
             v[3] = 1.f;
             
             vec4 wsv = gs2ws_transform * v;
             
             //cerr << "p = " << p << endl;
             //cerr << "v = " << v << endl;
             //cerr << "wsv = " << wsv << endl;
             
             const vec3 diff = vec3(wsv[0] - p[0], wsv[1] - p[1], wsv[2] - p[2]);
             const float r2 = dot(diff, diff);
             
             //cerr << "diff = " << diff << endl;
             //cerr << "length(diff) = " << r << endl;
             const float val = rdf.evaluate(r2);
             //cerr << "val = " << val << endl;
             
             float& res = rdf_volume->get_data(i,j,k);
             
             #pragma omp atomic
             res += val;
          }
        }
      }
    }
    
    clock.stop();
    cerr << "create_rdf_volume took " << clock.msecs() << " msecs." << endl;
    
    rdf_volume->ws2gs_scale = ws2gs_transform.diag3();
    rdf_volume->ws2gs_translate = ws2gs_transform.get_translate_vec3();
    rdf_volume->xyz_filebase = filebase;
    
    rdf_volume->compute_minmax();
    //rdf_volume->write();
  }
#endif

#if 1
  //AARONBAD -- this is for reading a special .dat file from Priya Vashishta's group
  //NOTE: the idea is to read a .raw/.nhdr first (this contains translate, scale info) and then read the .dat
  
  bool read_dat(string _filename = "") //, int flip_xz, float br = 1.f, float sr = 1.f, float sd = 0.f)
  {    
    if (_filename == "")
      _filename = filebase + ".dat";
    
    FILE* file = fopen(_filename.c_str(), "r");
    
    if (!file)
    {
      cerr << "Failed to load " << _filename << endl;
      return false;
    }
    
    int npts;

    float r0, r1;
    fscanf(file,"%f %f\n",&r0, &r1);

    fscanf(file,"%d\n",&npts);
    cerr << "Reading .dat file with " << npts << " atoms." << endl;
    
    int num_gs_atoms = npts;

    int* pid_of_idx = new int[num_gs_atoms];
    bond_adjacencies = new vector<int>[num_gs_atoms];

    float x,y,z,w;
	
		sort(special_atoms.begin(), special_atoms.end());
		cerr << "special_atoms.size() = " << special_atoms.size() << endl;
    
    int identifier;
    int blah;
		int sai=0;
    int points_read = 0;
    ws_atoms.resize(num_gs_atoms);

    vec4 min_atom = vec4(std::numeric_limits<float>::max());
    vec4 max_atom = vec4(std::numeric_limits<float>::min());

    for(int i = 0; i < num_gs_atoms; i++)
    {
			fscanf(file,"%d %d %f %f %f %d\n", &pid_of_idx[i], &identifier, &x,&y,&z, &blah);
			if (sai < special_atoms.size() && i == special_atoms[sai])
			{
				sai++;
				identifier = 0;
			}

      //if (pid_of_idx[i] != pid_of_idx[i-1] + 1)
      //  printf("Skipped something before %d\n", i);
      
      //NOTE: this is special-cased for Priya Vashishta's data. 
      //We might want to split this out into a separate file, instead of keeping in cubeutils.hxx
      {    
        if (identifier == 1)
          w = 8.f;  //oxygen
        else
        {
          w = 14.f; //silicon
        }
                         
        //cout << x << " " << y << " " << z << " " << w << endl;

        
        ws_atoms[i] = vec4f(x,y,z,w);
        //Matrix44::transform(gs_atoms[i], ws2gs_transform, ws_atoms[i]);
        					
				//if (!noshift)
				//{
					//cerr << "Shifting gs_atoms by half a voxel." << endl;
					//gs_atoms[i] += vec4(.5f, .5f, .5f, 0.f);
				//}

        //cout << "vertex = " << gs_atoms[i] << endl;
        
        min_atom = min(ws_atoms[i], min_atom);
        max_atom = max(ws_atoms[i], max_atom);
        
        points_read++;
      }
    }

    //NOTE, even if we don't have volume data, we still want a "grid" scale from which to make macrocells
    //this is determined by command-line options

	  //if ws2gs_scale != 0, then user specifies a transform and dims. Else, we need to create one
    if (ws2gs_scale[0] == 0.f)
    {
      cerr << "Warning, read_xyz() was called with no volume or transform specified. Generating one." << endl;
      
      vec4 pad_vector = (max_atom - min_atom) * (rdf_pad_percent * .01f * .5f);
      cerr << "min_atom = " << min_atom << endl;
      cerr << "max_atom = " << max_atom << endl;
      cerr << "rdf_pad_percent = " << rdf_pad_percent << endl;
      vec4 max_atom_pad = max_atom + pad_vector;
      vec4 min_atom_pad = min_atom - pad_vector;
      
      vec4 pad_extents = max_atom_pad - min_atom_pad;
      cerr << "pad_extents = " << pad_extents << endl;
      
      if (volume_dims[0])
      {
        voxels_per_angstrom = volume_dims[0] / pad_extents[0];
        ws2gs_scale[0] = voxels_per_angstrom;
      }

      cerr << "voxels_per_angstrom = " << voxels_per_angstrom << endl;
      //if volume_dims == 0, we need to determine them
      volume_dims[0] = int(ceilf(pad_extents[0] * voxels_per_angstrom));
      volume_dims[1] = int(ceilf(pad_extents[1] * voxels_per_angstrom));
      volume_dims[2] = int(ceilf(pad_extents[2] * voxels_per_angstrom));
      cerr << "volume_dims = " << volume_dims << endl;
      
      ws2gs_scale = vec4(voxels_per_angstrom);
      ws2gs_scale[3] = 1.f;
      
      ws2gs_translate = -min_atom_pad * ws2gs_scale;
      ws2gs_translate[3] = 0.f;
      
    }
    else //ws2gs_scale and ws2gs_translate have been specified
    {
    
      if (ws2gs_scale[1] == 0.f)
			  ws2gs_scale[1] = ws2gs_scale[0];
		  if (ws2gs_scale[2] == 0.f)
			  ws2gs_scale[2] = ws2gs_scale[0];
			  
			voxels_per_angstrom = ws2gs_scale[0];
    }
			
    cerr << "ws2gs_translate = " << ws2gs_translate << endl;
    cerr << "ws2gs_scale = " << ws2gs_scale << endl;
    
    ws2gs_transform.identity();
    ws2gs_transform.m[0][0] = ws2gs_scale[0];
    ws2gs_transform.m[1][1] = ws2gs_scale[1];
    ws2gs_transform.m[2][2] = ws2gs_scale[2];
    ws2gs_transform.m[0][3] = ws2gs_translate[0];
    ws2gs_transform.m[1][3] = ws2gs_translate[1];
    ws2gs_transform.m[2][3] = ws2gs_translate[2];
    ws2gs_transform.scrub();
    
    cerr << "ws2gs_transform = " << ws2gs_transform << endl;
    
    Matrix44::inverse(gs2ws_transform, ws2gs_transform);

    cerr << "gs2ws_transform = " << gs2ws_transform << endl;
   
    cerr << "Creating gs_atoms..." << endl;
    int sqrt_size = int(ceilf(sqrt(ws_atoms.size())));
#ifdef BAS_TEXTURE_BUFFER
    gs_atoms.reserve(ws_atoms.size());
#else
    gs_atoms.reserve(sqrt_size * sqrt_size);
#endif
    gs_atoms.resize(ws_atoms.size());
          
    for(int i=0; i<ws_atoms.size(); i++)
    {
      gs_atoms[i] = transform_atom(ws2gs_transform, ws_atoms[i]);
      gs_atoms[i] += vec4(.5f, .5f, .5f, 0.f);
    }
    
/*
    num_gs_atoms = points_read;
    
    cerr << "num_atoms = " << num_gs_atoms << endl;
    cerr << "done!" << endl;   
    
    build_gs_atoms();
    */

    //this will build bond_adjacencies
    if (cus->stick_distance > 0)
    { 
      const float stick_distance_squared = cus->stick_distance * cus->stick_distance;

      //create reverse mapping for pid_to_idx
      const int num_pids = pid_of_idx[gs_atoms.size()-1] + 1;
      cerr << "num_pids = " << num_pids << endl;
      int* idx_of_pid = new int[num_pids];
      memset(idx_of_pid, 0, num_pids * sizeof(int));
      for(int i=0; i<gs_atoms.size(); i++)
        idx_of_pid[pid_of_idx[i]] = i;

      //read in the sticks
      cerr << "reading in sticks from .dat file..." << endl;
      int i=0;
      int pid, stick_pid;
      char rest[32];

      const int BSIZE = 1024;
      char buffer[BSIZE];

      for(int i=0; i<gs_atoms.size(); i++)
      {
        fgets(buffer, BSIZE, file);

        sscanf(buffer, " %d", &pid);

        //if (pid != i+1)
        //  printf("\n\npid = %d, i=%d\n", pid, i);

        //printf("buffer = %s\n", buffer);

        if (idx_of_pid[pid] != i)
          printf("Error! idx_of_pid[%d] = %d, but we expected %d\n", pid, idx_of_pid[pid], i);

        vector<int>& idxsticks = bond_adjacencies[i];
        vec4 pv1 = gs_atoms[i];

        string sval = "";
        bool we_care = false;
        char last_c = ' ';
        for(char* c = buffer; ; c++)
        {
          if (*c == ' ' || *c == '\n')
          { 
            if (last_c != ' ')
            {
              int val = atoi(sval.c_str());
              if (we_care && val > 0)
              {
                vec4 pv2 = gs_atoms[idx_of_pid[val]];
  						  vec4f diff = pv2 - pv1;
							  diff *= diff;
							  float dist = diff[0] + diff[1] + diff[2];

                if (dist < stick_distance_squared)
                {
                  //printf(" %d ", val);
                  idxsticks.push_back( idx_of_pid[val] );
                }
              }
              we_care = !we_care;
              
              sval = "";

              if (*c == '\n')
                break;
            }
          }
          else if (*c != ' ')
            sval += *c;

          last_c = *c;
        }
        //printf("\n");
      }
    }

    cerr << "done!" << endl;
   
    fclose(file);
    
    return true;
  }
#endif

  void create_uchar_volume()
  {
    Array3<float>* current_volume = (Array3<float>*)volumes[cvi];
    
    Array3<unsigned char>* uchar_volume = new Array3<unsigned char>;
    volumes.push_back(uchar_volume);

    Clock clock;
    clock.start();
  
    const int num_voxels = current_volume->dims[0] * current_volume->dims[1] * current_volume->dims[2];
    uchar_volume->resize(current_volume->dims[0], current_volume->dims[1], current_volume->dims[2]);
    uchar_volume->type = "uchar";
    //uchar_volume->filebase = filebase + "_uchar";
    uchar_volume->filebase = filebase;
    uchar_volume->xyz_filebase = current_volume->xyz_filebase;
    float vrscale = 255.f / (current_volume->range[1] - current_volume->range[0]);
#pragma omp parallel for
    for(int i=0; i<num_voxels; i++)
      uchar_volume->data[i] = (unsigned char)((current_volume->data[i] - current_volume->range[0]) * vrscale);

    clock.stop();
    cerr << "create_uchar_volume took " << clock.msecs() << " msecs." << endl;

		if (cus)
			cus->display_uchar_volume = true;
      
    cvi++;
  }
	
	template<int NBALLS>
	struct GPBalls
	{
		vec4f balls[NBALLS];
	};
	
	template<int NBALLS>
	void build_gridded_particles(Array3<BuildMacrocell>& build_macrocells)
	{
		
		cerr << "build_gridded_particles<" << int(NBALLS) << ">()" << endl; 
		
		Array3< GPBalls<NBALLS> > gp(mc_dims[0], mc_dims[1], mc_dims[2]);
		gp.allocate();
		memset(&gp.data[0], gp.get_data_size(), 0);
		
#pragma omp parallel for		
    for(int k=0; k<mc_dims[2]; k++)
      for(int j=0; j<mc_dims[1]; j++)
        for(int i=0; i<mc_dims[0]; i++)
        {
          BuildMacrocell& cc = build_macrocells.get_data(i,j,k);
					GPBalls<NBALLS>& gc = gp.get_data(i,j,k);
					
					for(int b=0; b<cc.balls.size(); b++)
						gc.balls[b] = cc.balls[b];					
				}
	
		cerr << "done!" << endl;
		
		gp.filebase = filebase + "_gp";
		gp.write_raw();
		
	}	

  void write_volume()
  {
    Array3Base* current_volume = volumes[cvi];
    if (ws2gs_transform.is_cartesian())
    {
      current_volume->ws2gs_scale = ws2gs_transform.diag3();
      current_volume->ws2gs_translate = ws2gs_transform.get_translate_vec3();
    }
    current_volume->write();
  }
  
  bool write_xyz()
  {
    FILE* fout = fopen(string(filebase + ".xyz").c_str(), "w");
    if (!fout) return false;
    
    if (ws_atoms.size() == 0)
    {
      cerr << "Fatal: ws_atoms.size() = 0 in write_xyz()" << endl;
      return false;
    }
    
    fprintf(fout, "%d\n", int(ws_atoms.size()));
    //fprintf(fout, "nanovol02\n");
    //fprintf(fout, "%s: scale %f %f %f translate %f %f %f\n", description.c_str(), gs2ws_transform.m[0][0], gs2ws_transform.m[1][1], gs2ws_transform.m[2][2], gs2ws_transform.m[3][0], gs2ws_transform.m[3][1], gs2ws_transform.m[3][2]);
    
    for(int i=0; i<ws_atoms.size(); i++)
    {
      const vec4f& atom = ws_atoms[i];
      fprintf(fout, "%s %f %f %f\n", atom_str[int(atom[3])], atom[0], atom[1], atom[2]);
    }
    
    fclose(fout);
  }

  void count_atom_types()
  {
    int count[128];
    
    for(int i=0; i<128; i++)
      count[i] = 0;
      
    for(int i=0; i<ws_atoms.size(); i++)
      count[int(ws_atoms[i][3])]++;
    
    for(int i=0; i<128; i++)
    {
      if (count[i] > 0)
        atom_type_counts.push_back(vec2(i, count[i]));
    }
  }
  
  Array3<vec2uc>* create_minmax_macrocells()
  {
    Clock clock;
    clock.start();
    
    Array3<vec2uc>* minmax_macrocells;
    if (get_current_volume()->minmax_macrocells)
    {
      minmax_macrocells = (Array3<vec2uc>*)get_current_volume()->minmax_macrocells;
      return minmax_macrocells;
    }
    else
    { 
      minmax_macrocells = new Array3<vec2uc>();
      get_current_volume()->minmax_macrocells = minmax_macrocells;
    }
    
    minmax_macrocells->resize(mc_dims[0], mc_dims[1], mc_dims[2]);
    
    const vec2& float_data_range = get_current_volume()->range;
    const float float_data_scale = float(255.0 / double(float_data_range[1] - float_data_range[0]));
    
    const vec3i mcdm1 = mc_dims - vec3i(1);

    #pragma omp parallel for
    for(int k=0; k<mc_dims[2]; k++)
      for(int j=0; j<mc_dims[1]; j++)
        for(int i=0; i<mc_dims[0]; i++)
        {
          vec2uc& mm = minmax_macrocells->get_data(i,j,k);
          const int stencil_lo = 0;
          const int stencil_up = 1;

          int istart = i * voxels_per_macrocell - stencil_lo;
          int jstart = j * voxels_per_macrocell - stencil_lo;
          int kstart = k * voxels_per_macrocell - stencil_lo;
          int iend = istart + voxels_per_macrocell + stencil_up;
          int jend = jstart + voxels_per_macrocell + stencil_up;
          int kend = kstart + voxels_per_macrocell + stencil_up;

          istart = max(istart, 0);
          jstart = max(jstart, 0);
          kstart = max(kstart, 0);
          iend = min(iend, mcdm1[0]);
          jend = min(jend, mcdm1[1]);
          kend = min(kend, mcdm1[2]);

          mm[0] = 255;
          mm[1] = 0;
          
          for(int kk=kstart; kk<=kend; kk++)
            for(int jj=jstart; jj<=jend; jj++)
              for(int ii=istart; ii<=iend; ii++)
              {
                float val = get_current_volume()->get_element(ii,jj,kk);
                const unsigned char ucval = (val - float_data_range[0]) * float_data_scale;
                mm[0] = min(mm[0], ucval);
                mm[1] = max(mm[1], ucval);
              }
        }
        
    clock.stop();
    cerr << "create_minmax_macrocells took " << clock.msecs() << " msecs." << endl;
    
    return minmax_macrocells;
  }
  
  void update_tf_macrocells(TransferFunction* tf)
  {
    Array3<vec2uc>* minmax_macrocells = create_minmax_macrocells();
  
    Clock clock;
    clock.start();

    #pragma omp parallel for
    for(int k=0; k<mc_dims[2]; k++)
      for(int j=0; j<mc_dims[1]; j++)
        for(int i=0; i<mc_dims[0]; i++)
        {
          const vec2uc& mm = minmax_macrocells->get_data(i,j,k);
          vec4f& mc = macrocells.get_data(i,j,k);
          mc[0] = IDXCAST(tf->rangeMap[mm[0]][mm[1]]);     //a value between 0-255
        }
        
    clock.stop();
    cerr_dbg( "update_tf_macrocells took " << clock.msecs() << " msecs." << endl);
  }
   
  
  void create_macrocells(bool build_gpu_data = true)
  {
    //cerr << "In create_macrocells()" << endl;
    
    //first, determine macrocell size.
    if (voxels_per_macrocell == 0)
      voxels_per_macrocell = 4;

    cerr_dbg(  "voxels_per_macrocell = " << voxels_per_macrocell << endl  ); 

    const float stick_distance_scale_squared = voxels_per_angstrom * voxels_per_angstrom;
    cerr_dbg( "stick_distance_scale_squared = " << stick_distance_scale_squared << endl );

    cerr_dbg(  "volume_dims = " << volume_dims << endl  ); 
    
    mc_dims = vec3i(int(volume_dims[0] / float(voxels_per_macrocell)), 
                    int(volume_dims[1] / float(voxels_per_macrocell)), 
                    int(volume_dims[2] / float(voxels_per_macrocell)));

    cerr << "macrocell mc_dims = " << mc_dims << endl; 
    
    for(int i=0; i<3; i++)
    {
      if (mc_dims[i] * voxels_per_macrocell < volume_dims[i]) 
        mc_dims[i]++;
    }
                
    const vec4f mcdims_minus_one = vec4f(mc_dims[0] - 1.f, mc_dims[1] - 1.f, mc_dims[2] - 1.f, 4096.f);
    cerr_dbg(  "mcdims_minus_one = " << mcdims_minus_one << endl  ); 
    
    Array3<BuildMacrocell> build_macrocells(mc_dims[0], mc_dims[1], mc_dims[2]);
   
    const vec4 inv_mc_mul = vec4(1.f / float(voxels_per_macrocell), 1.f / float(voxels_per_macrocell), 1.f / float(voxels_per_macrocell), 1.f);
    num_indices = 0;
    int ball_indices = 0;
    int stick_indices = 0;
    cerr_dbg(  "inv_mc_mul = " << inv_mc_mul << endl  ); 
    
    cerr << "building balls" << endl;
    
    Clock clock;
    clock.start();
    
    //balls
    for(int v=0; v<gs_atoms.size(); v++)
    {
      //vec4f vert = atoms[v];   //note, this is in world (angstrom) space
      vec4f gs_vert = gs_atoms[v];
      float radius = cus->atom_colors[int(gs_vert[3])][3] * voxels_per_angstrom;
      //if (cus->display_analytical_volume)
      //  radius *= 1.5f
      gs_vert[3] = 1.f;
      
      vec4 bnr = gs_vert * inv_mc_mul;
      //cerr << "bnr[" << v << "] = " << bnr << endl;
      build_macrocells.get_data(int(bnr[0]), int(bnr[1]), int(bnr[2])).balls_no_radius.push_back(v);
      
      vec4 ball_min = gs_vert - vec4(radius, radius, radius, 0.0);
      vec4 ball_max = gs_vert + vec4(radius, radius, radius, 0.0);
      
      ball_min *= inv_mc_mul;
      ball_max *= inv_mc_mul;
      
      ball_min = max(ball_min, vec4(0,0,0,0));
      ball_max = min(ball_max, mcdims_minus_one);
      
      for(int k=int(ball_min[2]); k<=int(ball_max[2]); k++)
        for(int j=int(ball_min[1]); j<=int(ball_max[1]); j++)
          for(int i=int(ball_min[0]); i<=int(ball_max[0]); i++)
          {
            //cerr << "i=" << i << ", j=" << j << ", k = " << k << endl; 
            BuildMacrocell& mc = build_macrocells.get_data(i,j,k);
            mc.balls.push_back(v);
            ball_indices++;
          }
    }
    
    clock.stop();
    cerr << "balls took " << clock.msecs() << " msecs." << endl;
    
    cerr_dbg(  "finished balls, balls_no_radius #indices = " << gs_atoms.size() << endl  ); 
    cerr_dbg(  "finished balls, #indices = " << ball_indices << endl  ); 
    
    cerr_dbg(  "stick_distance = " << cus->stick_distance << endl  ); 
    
    clock.start();
    
    //create the global sticks, i.e. bond_adjacencies
    if (!bond_adjacencies)    //if no bond adjacencies
    {
      cerr << "no bond adjacencies, found, building sticks." << endl;
      bond_adjacencies = new vector<int>[gs_atoms.size()];
      for(int a=0; a<ws_atoms.size(); a++)
        bond_adjacencies[a].reserve(8);

      const float mc_stick_distance = cus->stick_distance * voxels_per_angstrom / voxels_per_macrocell;  
      cerr << "mc_stick_distance = " << mc_stick_distance << endl;
      
      #pragma omp parallel for
	    for(int v0=0; v0<ws_atoms.size(); v0++)
	    {
        vec4 gs_vert = gs_atoms[v0];
        vec4 ms_vert = gs_vert * inv_mc_mul;

        
        vec4 gsv_min = ms_vert - vec4(mc_stick_distance);
        vec4 gsv_max = ms_vert + vec4(mc_stick_distance);
        
        gsv_min = max(gsv_min, vec4(0.f));
        gsv_max = min(gsv_max, mcdims_minus_one);
        
        //cerr << "gsv_min = " << gsv_min << endl;
        //cerr << "gsv_max = " << gsv_max << endl;
        
        for(int kk=int(gsv_min[2]); kk<=int(gsv_max[2]); kk++)          
          for(int jj=int(gsv_min[1]); jj<=int(gsv_max[1]); jj++)      
            for(int ii=int(gsv_min[0]); ii<=int(gsv_max[0]); ii++)  
            {
              BuildMacrocell& cc1 = build_macrocells.get_data(ii,jj,kk);

              vector<int>& bnr1 = cc1.balls_no_radius;
              
              for(int bb=0; bb<bnr1.size(); bb++)
	            {
                const int v1 = bnr1[bb];
                if (v0 >= v1) continue;
                
                const vec4& p0 = ws_atoms[v0];
                const vec4& p1 = ws_atoms[v1];
                
                vec4 diff = p0 - p1;
                const float dist2 = dot3(diff, diff);
                if (dist2 < cus->get_stick_distance_squared(int(p0[3]), int(p1[3])))
                {
                  #pragma omp critical
                  bond_adjacencies[v0].push_back(v1);
                }
              }
            }
          }
    }
    clock.stop();

    cerr << "bond_adjacencies took " << clock.msecs() << " msecs." << endl;
    
    clock.start();
   
    //sticks
    if (cus->build_sticks)
    {
      //if (bond_adjacencies)  //build ball & stick macrocells from global ball & stick indices, e.g. for Priya .dat files
      {
        #pragma omp parallel for
		    for(int v1=0; v1<gs_atoms.size(); v1++)
		    {
          for(int v1a=0; v1a<bond_adjacencies[v1].size(); v1a++)
          {
            int v2 = bond_adjacencies[v1][v1a];
            if (v2 <= v1) continue;
            
            vec4 pv1 = ws_atoms[v1];
				    vec4 pv2 = ws_atoms[v2];
      
				    vec4f diff = pv2 - pv1;
				    float dist2 = dot3(diff, diff);
            
            if (dist2 < cus->get_stick_distance_squared(int(pv1[3]), int(pv2[3])))
            {
              //pv1[3] = 1.f;
              //pv2[3] = 1.f;  
              pv1 = gs_atoms[v1];
              pv2 = gs_atoms[v2];            
				      vec4 pv1_min = min(pv1, pv2);
				      vec4 pv1_max = max(pv1, pv2);
		          pv1_min *= inv_mc_mul;
		          pv1_max *= inv_mc_mul;
				      pv1_min = max(pv1_min, vec4(0,0,0,0));
				      pv1_max = min(pv1_max, mcdims_minus_one);
              
              for(int k=int(pv1_min[2]); k<=int(pv1_max[2]); k++)
		            for(int j=int(pv1_min[1]); j<=int(pv1_max[1]); j++)
		              for(int i=int(pv1_min[0]); i<=int(pv1_max[0]); i++)
		              {
		                BuildMacrocell& cc = build_macrocells.get_data(i,j,k);
		                
		                #pragma omp critical
		                {
            				  cc.sticks.push_back(v1);
						          cc.sticks.push_back(v2);
						          stick_indices += 2;
						        }
                  }
            }
          }
        }
      }
    }
    
    clock.stop();
    cerr << "sticks took " << clock.msecs() << " msecs." << endl;
    cerr_dbg(  "finished sticks, indices = " << stick_indices << endl  );     
    
    num_indices = ball_indices + stick_indices;

#if 1
    //AARONBAD -- gridded particles for Sidharth
    if (write_gridded_particles)
    {
			//from build_macrocells, create a .nhdr file with each voxel == N particles
			//we can then convert this to .idx
	
	    clock.start();
	
			int maxBallSize = 0;
	    for(int k=0; k<mc_dims[2]; k++)
	      for(int j=0; j<mc_dims[1]; j++)
	        for(int i=0; i<mc_dims[0]; i++)
	        {
	          BuildMacrocell& cc = build_macrocells.get_data(i,j,k);
						maxBallSize = max<int>(maxBallSize, int(cc.balls.size()));
					}
			
			cerr << "maxBallSize = " << maxBallSize << endl;
			
			if (maxBallSize > 16)
				build_gridded_particles<32>(build_macrocells);
			else if (maxBallSize > 8)
					build_gridded_particles<16>(build_macrocells);
			else if (maxBallSize > 4)
					build_gridded_particles<8>(build_macrocells);
			else
					build_gridded_particles<4>(build_macrocells);
      
	    clock.stop();
	    cerr << "write_gridded_particles took " << clock.msecs() << " msecs." << endl;
			exit(0);
    }
#endif

    if (build_gpu_data)
    {
      clock.start();

      //AARONBAD -- this should really be put in a separate "build GPU macrocells" function.
      //the reason: we don't want to have to build min-max tree for some ball&stick analysis

      //build macrocells for GPU now.
      
      cerr_dbg(  "total num_indices = " << num_indices << endl  ); 
      int sqrt_num_indices = int(ceilf(sqrtf(num_indices)));

#ifdef BAS_TEXTURE_BUFFER
      indices = new indextype[num_indices];
#else
      indices = new indextype[sqrt_num_indices*sqrt_num_indices];
#endif
      
      macrocells.resize(mc_dims[0], mc_dims[1], mc_dims[2]);
      
      int cidx = 0;
      for(int k=0; k<mc_dims[2]; k++)
        for(int j=0; j<mc_dims[1]; j++)
          for(int i=0; i<mc_dims[0]; i++)
          {
            BuildMacrocell& cc = build_macrocells.get_data(i,j,k);
            vec4f& mc = macrocells.get_data(i,j,k);
            
            //AARONBAD -- this should contain the min-max. Build min-max tree here
            //mc[0] = IDXCAST(tf_macrocells.get_data(i,j,k));     //a value between 0-255
            
            //for grid:
            //M, balls, sticks, last_stick
            //without volume:
            //rbfs, balls, sticks, last_stick

            if (cc.balls.size() == 0 && cc.sticks.size() == 0)
            {
              int tmp = -1;
              mc[1] = IDXCAST(tmp);
            }
            else 
            {
              mc[1] = IDXCAST(cidx);
              for(int b=0; b<cc.balls.size(); b++)
              {
                indices[cidx++] = cc.balls[b];
              }
              //cerr_dbg(  "balls.size() = " << cc.balls.size() << endl  ); 
              
              mc[2] = IDXCAST(cidx);
              for(int s=0; s<cc.sticks.size(); s++)
              {
                indices[cidx++] = cc.sticks[s];
              }
              
              //cerr_dbg(  "sticks.size() = " << cc.sticks.size() << endl  ); 
              
              mc[3] = IDXCAST(cidx);            
            }
          }
          
      clock.stop();
      cerr << "gpu_data took " << clock.msecs() << " msecs." << endl;
    }
  }

  //Mij = cos(exp (Zi Zj / (ri - rj)^2))
  void compute_fd2(double* w_array, double* fd_array, float Nc, float a, float wmin, float wmax, int wbins, int i)
  {
    const int off = i * wbins;
    
    for(int i=0; i<wbins; i++)
      fd_array[wbins + off] = 0.f;
        
#pragma omp parallel for
    for(int wb=0; wb<wbins; wb++)
    {
      fd_array[wb + off] = 0.0;
      
      double w = wmin + wb * double(wmax - wmin) / double(wbins);
      
      w_array[wb] = w;
      
      for(int i=0; i<ws_atoms.size(); i++)
      {
        const vec4 vi = ws_atoms[i];
        double zi = vi[3];
        
        //special case for diagonal elements (i == j)
        fd_array[wb + off] += zi*zi;
        
        //off-diagonal elements (i != j)
        for(int j=0; j<i; j++)
        {
          const vec4 vj = ws_atoms[j];
          double zj = vj[3];
          
          vec4 dv = vi - vj;
          double d = sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
          d = max<double>(d, 1e-6);
          
          fd_array[wb + off] += 2.0 * zi * zj * cos(w/d);
        }        
      }
      
      fd_array[wb + off] *= Nc * exp(-(w*w)/(2.0*a));
    }
  }
  
  //from zi * exp(-a*(x-xi)^2)
  void compute_fd1(double* w_array, double* fd_array, float Nc, float a, float wmin, float wmax, int wbins, int i)
  {
    const int off = i * wbins;
    
    for(int i=0; i<wbins; i++)
      fd_array[wbins + off] = 0.f;
    
#pragma omp parallel for
    for(int wb=0; wb<wbins; wb++)
    {
      fd_array[wb + off] = 0.0;
      
      double w = wmin + wb * double(wmax - wmin) / double(wbins);
      
      w_array[wb] = w;
      
      for(int i=0; i<ws_atoms.size(); i++)
      {
        const vec4 vi = ws_atoms[i];
        double zi = vi[3];
        
        //special case for diagonal elements (i == j)
        fd_array[wb + off] += zi*zi;
        
        //off-diagonal elements (i != j)
        for(int j=0; j<i; j++)
        {
          const vec4 vj = ws_atoms[j];
          double zj = vj[3];
          
          vec4 dv = vi - vj;
          double d = sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
          
          fd_array[wb + off] += 2.0 * zi * zj * cos(w*d);
        }        
      }
      
      fd_array[wb + off] *= Nc * exp(-(w*w)/(2.0*a));
    }
  }

  //FD4 = Z1*Cos[(Z1*E^(-a (d - d11)^2) + Z2*E^(-a (d - d12)^2) + 
  //               Z3*E^(-a (d - d13)^2) + Z4*E^(-a (d - d14)^2))/Z1]
  //FD4(d) = zi * cos(sum_j zj*exp(-a (d-d_ij)^2))
  void compute_fd4(double* w_array, double* fd_array, float Nc, float a, float wmin, float wmax, int wbins, int i)
  {
    const int off = i * wbins;
    
    for(int i=0; i<wbins; i++)
      fd_array[wbins + off] = 0.f;
    
#pragma omp parallel for
    for(int wb=0; wb<wbins; wb++)
    {
      fd_array[wb + off] = 0.0;
      
      double w = wmin + wb * double(wmax - wmin) / double(wbins);
      
      w_array[wb] = w;
      
      double inner = 0.0;
      
      for(int i=0; i<ws_atoms.size(); i++)
      {
        const vec4 vi = ws_atoms[i];
        double zi = vi[3];
        
        //off-diagonal elements (i != j)
        for(int j=0; j<ws_atoms.size(); j++)
        {
          const vec4 vj = ws_atoms[j];
          double zj = vj[3];
          
          vec4 dv = vi - vj;
          double d = sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
        
          inner += zj * exp(-a * (w - d)*(w - d));
        }        
        
        //special case for diagonal elements (i == j)
        fd_array[wb + off] += zi * cos(inner);
      }
    }
  }

  
  //from zi * exp(-a*(x-xi)^2/Zi)
  void compute_fd3(double* w_array, double* fd_array, float Nc, float a, float wmin, float wmax, int wbins, int i)
  {
    const int off = i * wbins;
    
    for(int i=0; i<wbins; i++)
      fd_array[wbins + off] = 0.f;
    
#pragma omp parallel for
    for(int wb=0; wb<wbins; wb++)
    {
      fd_array[wb + off] = 0.0;
      
      double w = wmin + wb * double(wmax - wmin) / double(wbins);
      
      w_array[wb] = w;
      
      for(int i=0; i<ws_atoms.size(); i++)
      {
        const vec4 vi = ws_atoms[i];
        double zi = max<double>(vi[3], 1e-6);
        
        //special case for diagonal elements (i == j)
        fd_array[wb + off] += zi * exp(-w*w*zi / (2*a)) * pow(zi, 5.0);
        
        //off-diagonal elements (i != j)
        for(int j=0; j<i; j++)
        {
          const vec4 vj = ws_atoms[j];
          double zj = max<double>(vj[3], 1e-6);
          
          vec4 dv = vi - vj;
          double d = sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
          
          fd_array[wb + off] += 2.0 * exp(-(w*w)*(zi+zj)/(4.0*a)) * pow(zi * zj, 2.5) * cos(w*d);
        }        
      }
      
      fd_array[wb + off] = pow(fd_array[wb + off], 0.2);
      fd_array[wb + off] *= Nc;
    }    
  }
  
  void plot_descriptor(string outfile, double* w_array, double* fd_array, double a, int wbins, bool gnuplot)
  {
    //write outfile
    FILE* fout = fopen(outfile.c_str(), "w");
    if (!fout)
      cerr << "Error writing to file " << outfile << endl;
    else
    {
      for(int i=0; i<wbins; i++)
        fprintf(fout, "%f %f\n", w_array[i], fd_array[i]);
    }
    fclose(fout);
    
    if (gnuplot)
    {
      //create a .gnu    
      fout = fopen(string(outfile + ".gnu").c_str(), "w");
      if (!fout)
        cerr << "Error writing to file " << string(outfile + ".gnu") << endl;
      else
      {
        fprintf(fout, "set terminal postscript eps color\n");
        //fprintf(fout, "set title \"%s, %d atoms\"\n", outfile.c_str(), ws_atoms.size());
        fprintf(fout, "set output \"%s.gnu.eps\"\n", outfile.c_str());
        fprintf(fout, "set xrange[0:8]\n");
        fprintf(fout, "plot \"%s\" using 1:2 with lines title \"%s, a=%d, %d atoms\" \n", outfile.c_str(), outfile.c_str(), int(a), ws_atoms.size());
      }
      fclose(fout);
    }
  }
  
  
};


struct RadialDensityDistribution
{
  CubeData* cd;
  
  float cutoff;
  int bins;
  int max_npts;
  int num_angles;
  bool stochastic, bondlines, write;
    
  struct DensityData
  {
    float min;
    float bonds_min;
    float max;
    float bonds_max;
    float density;
    int weight;
  };
  
  DensityData** densityData;
  
  RadialDensityDistribution()
  {
    bins = 200;
    cutoff = 2.f;
    max_npts = 0;
    num_angles = 10;
    stochastic = true;
    bondlines = true;
    write = true;
  }
  ~RadialDensityDistribution()
  {
    for(int c=0; c<cd->atom_type_counts.size(); c++)
      delete densityData[c];
    delete[] densityData;
  }
  
  float map_value(int atom, float r)
  {
    int bin = int(r * bins / cutoff);
    return densityData[atom][bin].min;
  }
  
  vec2f map_interval(int atom, float r_lo, float r_hi)
  {
    int bin_lo = int(r_lo * bins / cutoff);
    int bin_hi = int(ceilf(r_hi * bins / cutoff));
    
    vec2f range;
    range[0] = std::numeric_limits<float>::max();
    range[1] = std::numeric_limits<float>::min();

    for(int i=bin_lo; i<=bin_hi; i++)
    {
      range[0] = min(range[0], densityData[atom][i].min);
      range[1] = max(range[1], densityData[atom][i].bonds_min);
    }
    
    return range;
  }
  

  //AARONBAD
  //this VERY POORLY NAMED function computes radial density distributions of rho or psi 
  int compute_density_distributions()
  {
    Clock clock;
    
    clock.start();
    
    if (cd->atom_type_counts.size() == 0)
      cd->count_atom_types();
    
    /*
    if (cd->atom_type_counts.size() == 0)
    {
      cd->count_atom_types();
      if (cd->atom_type_counts.size() == 0)
      {
        cerr_dbg(  "radial_distribution_function error: no atoms counted!" << endl  ); 
        return 0;
      }
    }
     */
    
    cerr_dbg(  "computing density distributions..." << endl  ); 
    
    cerr_dbg(  "atom_type_counts.size() = " << cd->atom_type_counts.size() << endl  ); 
    
    densityData = new DensityData*[cd->atom_type_counts.size()];
    
    Array3<float>* float_volume = (Array3<float>*)cd->volumes[0];
    
    printf("float_volume range is %f, %f\n", float_volume->range[0], float_volume->range[1]);
    
    for(int c=0; c<cd->atom_type_counts.size(); c++)
    {
      densityData[c] = new DensityData[bins];
      
      for(int i=0; i<bins; i++)
      {
        densityData[c][i].density = 0.f;
        densityData[c][i].weight = 0;
        densityData[c][i].min = float_volume->range[1];
        densityData[c][i].bonds_min = float_volume->range[1];
        densityData[c][i].max = float_volume->range[0];
        densityData[c][i].bonds_max = float_volume->range[0];
      }
    }
    
    char atommap[128];
    for(int c=0; c<cd->atom_type_counts.size(); c++)
    {
      atommap[ cd->atom_type_counts[c][0] ] = c;
    }
    
    vec4 vcutoff = vec4(cutoff, cutoff, cutoff, cutoff);
    
    //now we have a float volume, and a .xyz file
    
    if (max_npts == 0)
      max_npts = cd->ws_atoms.size();
      
    if (cd->ws_atoms.size() < max_npts)
      max_npts = cd->ws_atoms.size();
    
    vec4 scale = cd->ws2gs_transform * vec4(1.f,0.f,0.f,0.f);
    scale[3] = 0.f;
    const float ws2gs_radius_scale = length(scale);
    cerr_dbg(  "ws2gs_radius_scale = " << ws2gs_radius_scale << endl  ); 
    
    cerr_dbg(  "max_npts = " << max_npts << endl  ); 
    cerr_dbg(  "num_angles = " << num_angles << endl  ); 
    
    if (stochastic)
    {
      cerr_dbg(  "stochastic line sampling..." << endl  ); 
      clock.start();

      for(int i=0; i<max_npts; i++)
      {
        vec4 atom = cd->ws_atoms[i];
       
        cerr_dbg(  "atom = " << atom << endl  ); 
        
        const int am = atommap[int(atom[3])];
        
        atom[3] = 1.f;    
        vec4 gs_atom = cd->ws2gs_transform * atom;

        cerr_dbg(  "atom " << i << ", gs_atom = " << gs_atom << endl  ); 

        
        DensityData* dd = densityData[am];
        
        for(int j=0; j<num_angles; j++)
        {
          //pick a random theta, phi (polar coordinates) between 0 and 2pi
          const double _2pi = 3.141592653589793f * 2.f;
          const double theta = drand48() * _2pi;
          const double phi = drand48() * _2pi;
          
          //find the vector corresponding to this
          vec4 v;
          v[0] = sin(theta) * cos(phi);
          v[1] = sin(theta) * sin(phi);
          v[2] = cos(phi);
          v[3] = 0.f;
          
          //normalize(v);
          
          v[3] = 1.f;
          
          v *= ws2gs_radius_scale * cutoff / bins;
          vec4 gs_sample = gs_atom; // - vec4(.5f, .5f, .5f, 0.f);
          
          //for each r from 0 to bins, 
          for(int k=0; k<bins; k++)
          {
            //sample the volume at radial distance k
            gs_sample += v;
            
            //vec3 gspos = vec3(gs_sample[0], gs_sample[1], gs_sample[2]);
            float density = float_volume->trilerp(((vec3&)gs_sample));
            
            //cerr_dbg(  "k = " << k << ", density = " << density << endl  ); 
            
            dd[k].density += density;
            dd[k].weight++;
            dd[k].min = min(dd[k].min, density);
            dd[k].max = max(dd[k].max, density);
          }
        }
        
        clock.stop();
        cerr_dbg(  "stochastic line sampling done in " << clock.msecs() << " ms." << endl  ); 
        
      }
    }
    
    if (bondlines)
    {
      clock.start();
      
      cerr_dbg(  "bond line sampling..." << endl  ); 
      
      //like stochastic (interpolates), but samples along bond lines
      
      cd->create_macrocells();
      
      for(int i=0; i<max_npts; i++)
      {
        vec4 atom = cd->ws_atoms[i];
        
        const int am = atommap[int(atom[3])];
        
        atom[3] = 1.f;    
        vec4 gs_atom = cd->ws2gs_transform * atom;
        
        //cerr_dbg(  "atom " << i << ", gs_atom = " << gs_atom;
        
        DensityData* dd = densityData[am];
        
        int num_lines = cd->bond_adjacencies[i].size();
        if (num_angles > 0)
          num_lines = min(num_lines, num_angles);
        
        //cerr_dbg(  " has " << num_lines << " adjacencies." << endl  ); 
        
        for(int j=0; j<num_lines; j++)
        {
          
          //if (j <= i) continue;
          
          
          //find the vector corresponding to this
          vec4 atom2 = cd->ws_atoms[ cd->bond_adjacencies[i][j] ];
          atom2[3] = 1.f;
          
          
          vec4 gs_atom2 = cd->ws2gs_transform * atom2;
          vec4 v = normalize(gs_atom2 - gs_atom);
          v[3] = 1.f;
          
          //cerr_dbg(  "adjacency " << j << ", gs_atom2 = " << gs_atom2 << endl  ); 
          
          //cerr_dbg(  "v = " << v << endl << endl  ); 
          
          v *= ws2gs_radius_scale * cutoff / bins;
          vec4 gs_sample = gs_atom; // - vec4(.5f, .5f, .5f, 0.f);
          
          //for each r from 0 to bins, 
          for(int k=0; k<bins; k++)
          {
            //sample the volume at radial distance k
            gs_sample += v;
            
            //vec3 gspos = vec3(gs_sample[0], gs_sample[1], gs_sample[2]);
            float density = float_volume->trilerp(((vec3&)gs_sample));
            
            //cerr_dbg(  "k = " << k << ", density = " << density << endl  ); 
            
            dd[k].density += density;
            dd[k].weight++;
            dd[k].bonds_min = min(dd[k].bonds_min, density);
            dd[k].bonds_max = max(dd[k].bonds_max, density);
          }
          
        }
        
        clock.stop();
        cerr_dbg(  "bond line sampling done in " << clock.msecs() << " ms." << endl  ); 
      }
      
      //normalize density over weight
      for(int c=0; c<cd->atom_type_counts.size(); c++)
        for(int i=0; i<bins; i++)
        {
          if (densityData[c][i].weight > 0)
            densityData[c][i].density /= densityData[c][i].weight;
          
          if (densityData[c][i].min == float_volume->range[1])
            densityData[c][i].min = 0.f;
          
          if (densityData[c][i].min > densityData[c][i].bonds_min)
            densityData[c][i].min = densityData[c][i].bonds_min;
        }
      
      cerr_dbg(  "done!" << endl  ); 
      clock.stop();
      cerr_dbg(  "Computing RDF's took " << clock.msecs() << " ms." << endl  ); 
      
      if (write)
      {        
        //write out text files
        for(int c=0; c<cd->atom_type_counts.size(); c++)
        {
          char pdffilename[512];
          sprintf(pdffilename, "%s_bulkdist_%s.txt", cd->filebase.c_str(), atom_str[cd->atom_type_counts[c][0]]);
          FILE* fout = fopen(pdffilename, "w");
          for(int i=0; i<bins; i++)
            fprintf(fout, "%d %f %f %f %f %f\n", i, i * cutoff / float(bins), densityData[c][i].density, densityData[c][i].min, densityData[c][i].bonds_min, densityData[c][i].max);
          fclose(fout);
          
          cerr_dbg(  "Wrote " << pdffilename << endl  ); 
        }
      }
      

    }
    
    return 1;
  
  }  
  
};

struct CubeDataSequence
{
  CubeUtilsSettings* default_cus;
  
  string filebase;
  string name;
  string description;

  vec3 domain_max;
  
  int current_step;   //current step in CubeDataSequence

  CubeDataSequence()
  {
    domain_max = vec3(0,0,0);
    default_cus = 0;

  }
  ~CubeDataSequence(){}

  void read_filename(string filenameString, int argc, char* argv[])
  { 
    string suffix;
    int pos = filenameString.find_last_of(".");
    if (pos > 0)
    {
      string suffix = filenameString.substr(pos);
      if (suffix == ".seq" || suffix == ".cds")
      {
        read_sequence_filename(filenameString, argc, argv);
        return;
      }
    }
    else
    {
      cerr << "no suffix found..." << endl;
    }

    CubeUtilsSettings* cus = new CubeUtilsSettings(default_cus);
    CubeData* cd = new CubeData();
    cd->cus = cus;    
    cd->parse_args(argc, argv);
    cus->parse_args(argc, argv);
    
    if (!cd->read_filename(filenameString))
    {
      cerr << "Error reading " << filenameString << endl;
      delete cd;
      delete cus;
    }
    else
    {
      //AARONBAD -- why do we do this?
      //cd->create_uchar_volume();
      //cd->write_uchar_volume();
      
      //for writing data after reads (e.g. RDF volumes, conversion)
      cd->parse_args_after_read(argc, argv);
      
      cd->create_macrocells();
      
      domain_max = max(domain_max, vec3(cd->volume_dims[0], cd->volume_dims[1], cd->volume_dims[2]));
      
      steps.push_back(cd);
    }

  }

  void read_sequence_filename(string filename, int argc, char* argv[])
  {
    string path;
    int slash = filename.find_last_of("/");
    if (slash != filename.npos)
      path = filename.substr(0, slash) + "/";
    else
      path = "";
    
    cerr << "Reading CubeDataSequence file " << filename << endl;
    cerr << "path = " << path << endl;
    
    char fbuffer[1024];
    FILE* fin = fopen(filename.c_str(), "r+");
    while(fscanf(fin, "%s\n", fbuffer) != EOF)
    {
      string seqfile = path + string(fbuffer);
      cerr << "in sequence file " << filename << ", component = " << seqfile << endl;
      read_filename(seqfile, argc, argv);
    }
    fclose(fin);
  }
  
  void plot_descriptor(string outfile, double* w_array, double* fd_array, double a, int wbins, bool gnuplot)
  {
    //write outfile
    FILE* fout = fopen(outfile.c_str(), "w");
    if (!fout)
      cerr << "Error writing to file " << outfile << endl;
    else
    {
      for(int s=0; s<steps.size(); s++)
      {
        fprintf(fout, "\n");
        for(int i=0; i<wbins; i++)
          fprintf(fout, "%f %f %f\n", float(s), w_array[i], (fd_array[s*wbins + i]));
      }
    }
    fclose(fout);
    
    if (gnuplot)
    {
      //create a .gnu    
      fout = fopen(string(outfile + ".gnu").c_str(), "w");
      if (!fout)
        cerr << "Error writing to file " << string(outfile + ".gnu") << endl;
      else
      {
        fprintf(fout, "set terminal postscript eps color\n");
        //fprintf(fout, "set title \"%s, %d atoms\"\n", outfile.c_str(), ws_atoms.size());      
        fprintf(fout, "set output \"%s.gnu.eps\"\n", outfile.c_str());
        //fprintf(fout, "set pm3d\n");
        fprintf(fout, "unset surface\n");
        fprintf(fout, "set samples 50, 50; set isosamples 50, 50\n");
        fprintf(fout, "set style line 100 lt 5 lw 0.25\n");
        fprintf(fout, "set pm3d hidden3d 100\n");
        fprintf(fout, "set view 50,135\n");
        fprintf(fout, "set pm3d implicit at s\n");
        fprintf(fout, "set pm3d scansbackward\n");
        fprintf(fout, "set xrange [0:%d]\n", steps.size());
        fprintf(fout, "set yrange [0:8]\n");
        //fprintf(fout, "set zrange [0:30]\n", steps.size());
        //fprintf(fout, "set yrange [%f:%f]\n", amplitude_min, amplitude_max);
        fprintf(fout, "splot \"%s\" using 1:2:3\n", outfile.c_str());
        //fprintf(fout, "replot\n", steps.size());
      }
      fclose(fout);
    }
  }
    
  

  
  vector<CubeData*> steps;
};


#endif
