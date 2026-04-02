#ifndef PARAMETERS_HH
#define PARAMETERS_HH


// ***** DON'T TOUCH *****
#define GAUSSIAN 0
// ***********************

// Initial data
#define INITIAL_DATA GAUSSIAN


// Diffusion constant
constexpr inline double ALPHA = 1.0;

// Adiabatic exponent
constexpr inline double GAMMA = 1.0;


// Left and right time-independent function values (Dirichlet BCs)
constexpr inline double U0 = 0.0;
constexpr inline double UL = 0.0;


// Number of spatial grid points
constexpr inline size_t NX = 1000;

// Spatial bounds are [0, L - L/N]
constexpr inline double L = 1.0;


// Number of time steps
constexpr inline size_t NT = 10000;

// Time step
constexpr inline double DT = 1.0e-08;


// Initial function center
constexpr inline double X0 = 0.3*L;

// Initial wavefunction half-width
constexpr inline double SIGMA = 0.03*L;


// Newton-Raphson tolerance to linear system solution (should be "small")
constexpr inline double NRTOL = 1.0e-12;

// Maximum number of Newton-Raphson iterations
constexpr inline size_t NRMAXIT = 100;


// Output frequency
constexpr inline size_t OUT_EVERY = 100;

// Output filename
#define FILENAME "Data.h5"

// Verbosity
constexpr inline bool VERBOSE = false;


#endif  // PARAMETERS_HH
