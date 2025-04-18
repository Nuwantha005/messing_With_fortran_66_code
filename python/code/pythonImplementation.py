import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
import math

class TurbomachinerySolver:
    """
    Solves turbomachinery flow using a streamline curvature method.

    This class encapsulates the data and methods translated from the
    legacy Fortran IV code. It reads input, performs iterative
    calculations involving spline fitting and interpolation, and
    outputs results.
    """
    def __init__(self):
        """Initializes solver parameters and data arrays."""
        # Constants
        self.RHO = 0.0  # Density (placeholder, read from input)
        self.GAM = 0.0  # Ratio of specific heats
        self.AR = 0.0   # Gas constant
        self.TEMP = 0.0 # Inlet temperature
        self.W = 0.0    # Rotational speed (rad/s)
        self.WT = 0.0   # Target weight flow
        self.ALM = 0.0  # Inlet swirl parameter
        self.TOLER = 0.0 # Convergence tolerance for streamline position
        self.PLOSS = 0.0 # Pressure loss parameter
        self.WTOLER = 0.0 # Convergence tolerance for weight flow
        self.MTHTA = 0  # Number of points for blade thickness definition
        self.NPRT = 0   # Print frequency for streamlines
        self.ITER = 0   # Max iterations / initial flag
        self.SFACT = 0.0 # Splitter factor
        self.ZSPLIT = 0.0 # Z-coordinate for splitter effect
        self.BETIN = 0.0 # Inlet blade angle (radians)
        self.RB = 0.0   # Base radius for blade angle calculation
        self.CORFAC = 0.0 # Correction factor for streamline update
        self.XN = 0.0   # Number of blades

        # Grid dimensions
        self.MX = 0     # Number of points along meridional direction
        self.KMX = 0    # Number of streamlines
        self.MR = 0     # Number of radial points for thickness input
        self.MZ = 0     # Number of axial points for thickness input

        # Input geometry and blade data
        self.ZS = None  # Shroud Z coordinates (MX)
        self.ZH = None  # Hub Z coordinates (MX)
        self.RS = None  # Shroud R coordinates (MX)
        self.RH = None  # Hub R coordinates (MX)
        self.THTA = None # Blade angle data (MTHTA)
        self.XT = None  # Blade thickness data (MTHTA)
        self.TN = None  # Blade thickness grid (MZ, MR)
        self.XZ = None  # Axial coordinates for thickness grid (MZ)
        self.XR = None  # Radial coordinates for thickness grid (MR)

        # Calculated variables (Arrays typically MX x KMX)
        self.AL = None  # Meridional angle
        self.BETA = None # Relative flow angle
        self.CAL = None # Cos(AL)
        self.CBETA = None # Cos(BETA)
        self.CURV = None # Meridional curvature
        self.DN = None  # Normal distance along quasi-orthogonals
        self.PRS = None # Static pressure
        self.R = None   # Radial coordinate of streamlines
        self.Z = None   # Axial coordinate of streamlines
        self.SM = None  # Meridional distance along streamlines
        self.SA = None  # Coefficient A in momentum equation
        self.SB = None  # Coefficient B in momentum equation
        self.SC = None  # Coefficient C in momentum equation
        self.SD = None  # Coefficient D in momentum equation
        self.SAL = None # Sin(AL)
        self.SBETA = None # Sin(BETA)
        self.TN_interp = None # Interpolated thickness (placeholder)
        self.TT = None  # Tangential thickness
        self.WA = None  # Absolute velocity magnitude (or axial component?) - Check usage
        self.WTR = None # Relative velocity magnitude on blade surface

        # Intermediate calculation arrays (typically size MX or KMX)
        self.AB = None
        self.AC = None
        self.AD = None
        self.BA = None
        self.DELBTA = None
        self.DRDM = None
        self.DTDR = None
        self.DTDZ = None
        self.DWMDM = None
        self.DWTDM = None
        self.WTFL = None

        # Control flags
        self.TYPE = 0   # Input type flag
        self.BCDP = 0   # BCD dump flag
        self.SRW = 0    # Unused in Python? (Appears related to COMMON)

        # Calculated constants
        self.CI = 0.0   # Inlet speed of sound
        self.CP = 0.0   # Specific heat at constant pressure
        self.EXPON = 0.0 # Exponent (1 / (GAM - 1.0))
        self.KMXM1 = 0

        # Iteration tracking
        self.runo = 0
        self.itno = 1
        self.error = 100000.0
        self.errori = 100000.0

        # Constants
        self.PI = math.pi
        self.ROOT2 = math.sqrt(2.0)
        self.DEG_TO_RAD = self.PI / 180.0
        self.RAD_TO_DEG = 180.0 / self.PI
        self.PSI_TO_PSF = 144.0
        self.IN_TO_FT = 1.0 / 12.0

    def _read_floats(self, n):
        """Helper to read n floats from a line of input."""
        return list(map(float, input().split()[:n]))

    def _read_ints(self, n):
        """Helper to read n integers from a line of input."""
        return list(map(int, input().split()[:n]))

    def read_input(self):
        """Reads input data corresponding to Fortran READ statements."""
        print("Reading input...")
        self.runo += 1

        # READ(5,1010) MX, KMX, MR, MZ, W, WT, XN, GAM, AR
        line1 = input().split()
        self.MX = int(line1[0])
        self.KMX = int(line1[1])
        self.MR = int(line1[2])
        self.MZ = int(line1[3])
        self.W = float(line1[4])
        self.WT = float(line1[5])
        self.XN = float(line1[6])
        self.GAM = float(line1[7])
        self.AR = float(line1[8])
        self.KMXM1 = self.KMX - 1

        print(f"Run No. {self.runo} Input Data Card Listing")
        print(f"{self.MX:5d}{self.KMX:5d}{self.MR:5d}{self.MZ:5d}"
              f"{self.W:10.4f}{self.WT:10.4f}{self.XN:10.4f}"
              f"{self.GAM:10.4f}{self.AR:10.4f}")

        # READ(5,1010) TYPE, BCDP, SRW, NULL, TEMP, ALM, RHO, TOLER, PLOSS, WTOLER
        line2 = input().split()
        self.TYPE = int(line2[0])
        self.BCDP = int(line2[1])
        self.SRW = int(line2[2])
        # NULL = int(line2[3]) # Ignored
        self.TEMP = float(line2[4])
        self.ALM = float(line2[5])
        self.RHO = float(line2[6])
        self.TOLER = float(line2[7])
        self.PLOSS = float(line2[8])
        self.WTOLER = float(line2[9])
        print(f"{self.TYPE:5d}{self.BCDP:5d}{self.SRW:5d}{0:5d}" # Assuming NULL was 0
              f"{self.TEMP:10.4f}{self.ALM:10.4f}{self.RHO:10.4f}"
              f"{self.TOLER:10.4f}{self.PLOSS:10.4f}{self.WTOLER:10.4f}")

        # READ(5,1010) MTHTA, NPRT, ITER, NULL, SFACT, ZSPLIT, BETIN, RB, CORFAC
        line3 = input().split()
        self.MTHTA = int(line3[0])
        self.NPRT = int(line3[1])
        self.ITER = int(line3[2])
        # NULL = int(line3[3]) # Ignored
        self.SFACT = float(line3[4])
        self.ZSPLIT = float(line3[5])
        self.BETIN = float(line3[6]) # Degrees initially
        self.RB = float(line3[7])
        self.CORFAC = float(line3[8])
        print(f"{self.MTHTA:5d}{self.NPRT:5d}{self.ITER:5d}{0:5d}" # Assuming NULL was 0
              f"{self.SFACT:10.4f}{self.ZSPLIT:10.4f}{self.BETIN:10.4f}"
              f"{self.RB:10.4f}{self.CORFAC:10.4f}")

        # Allocate arrays based on dimensions
        self._allocate_arrays()

        # Read geometry arrays
        self.ZS = np.array(self._read_floats(self.MX))
        print(*(f"{x:10.4f}" for x in self.ZS))
        self.ZH = np.array(self._read_floats(self.MX))
        print(*(f"{x:10.4f}" for x in self.ZH))
        self.RS = np.array(self._read_floats(self.MX))
        print(*(f"{x:10.4f}" for x in self.RS))
        self.RH = np.array(self._read_floats(self.MX))
        print(*(f"{x:10.4f}" for x in self.RH))

        # Unit conversion: inches to feet
        self.ZS *= self.IN_TO_FT
        self.ZH *= self.IN_TO_FT
        self.RS *= self.IN_TO_FT
        self.RH *= self.IN_TO_FT

        # Initialize flow field arrays based on TYPE
        if self.TYPE == 0:
            # Calculate inlet velocity assuming uniform flow
            inlet_area = self.PI * (self.RS[0]**2 - self.RH[0]**2) # Approx area, Fortran used mean radius
            inlet_area_corrected = self.PI * (self.RS[0] + self.RH[0]) * (self.ZS[0] - self.ZH[0]) # Closer to Fortran? Check geometry assumptions
            wa_inlet = self.WT / (self.RHO * inlet_area_corrected) # Fortran used WA(1,1) = WT/RHO/(ZS(1)-ZH(1))/3.14/(RS(1)+RH(1))

            for i in range(self.MX):
                dn_kmx = math.sqrt((self.ZS[i] - self.ZH[i])**2 + (self.RS[i] - self.RH[i])**2)
                for k in range(self.KMX):
                    frac = k / self.KMXM1
                    self.DN[i, k] = frac * dn_kmx
                    self.WA[i, k] = wa_inlet # Initial guess
                    self.Z[i, k] = frac * (self.ZS[i] - self.ZH[i]) + self.ZH[i]
                    self.R[i, k] = frac * (self.RS[i] - self.RH[i]) + self.RH[i]
            self.DN[:, 0] = 0.0 # Hub streamline DN is 0

        elif self.TYPE == 1:
            print("Reading BCD data for DN, WA, Z, R...")
            # Assumes BCD data is read row by row (Fortran reads column-major)
            # This part needs careful adaptation based on how BCDUMP/BCREAD work
            # and how the data is formatted. Placeholder:
            print("BCREAD for DN:")
            self.DN = np.array([self._read_floats(self.KMX) for _ in range(self.MX)])
            print("BCREAD for WA:")
            self.WA = np.array([self._read_floats(self.KMX) for _ in range(self.MX)])
            print("BCREAD for Z:")
            self.Z = np.array([self._read_floats(self.KMX) for _ in range(self.MX)])
            print("BCREAD for R:")
            self.R = np.array([self._read_floats(self.KMX) for _ in range(self.MX)])
            # Apply unit conversions if BCD data is in inches
            # self.DN *= self.IN_TO_FT
            # self.Z *= self.IN_TO_FT
            # self.R *= self.IN_TO_FT
            print("Finished reading BCD data.")
        else:
            raise ValueError(f"Unsupported TYPE: {self.TYPE}")

        # Read blade data
        self.THTA = np.array(self._read_floats(self.MTHTA))
        print(*(f"{x:10.4f}" for x in self.THTA))
        self.XT = np.array(self._read_floats(self.MTHTA))
        print(*(f"{x:10.4f}" for x in self.XT))

        # Read thickness grid
        self.TN = np.zeros((self.MZ, self.MR))
        for k in range(self.MR): # Read row by row (Fortran reads column major for TN(I,K))
             row = self._read_floats(self.MZ)
             self.TN[:, k] = row # Store column-wise to match Fortran TN(I,K) where I=MZ, K=MR
             print(*(f"{x:10.4f}" for x in row))

        self.XZ = np.array(self._read_floats(self.MZ))
        print(*(f"{x:10.4f}" for x in self.XZ))
        self.XR = np.array(self._read_floats(self.MR))
        print(*(f"{x:10.4f}" for x in self.XR))

        # Unit conversions for blade/thickness data
        self.TN *= self.IN_TO_FT
        self.XR *= self.IN_TO_FT
        self.XZ *= self.IN_TO_FT
        self.XT *= self.IN_TO_FT # Thickness data

        # Initialize other parameters and arrays
        self.SM.fill(0.0)
        self.BA = np.linspace(0.0, self.WT, self.KMX) # Weight flow distribution guess

        # Convert remaining input parameters
        self.TOLER *= self.IN_TO_FT
        self.RB *= self.IN_TO_FT
        self.ZSPLIT *= self.IN_TO_FT
        self.PLOSS *= self.PSI_TO_PSF
        self.BETIN *= -self.DEG_TO_RAD # Convert degrees to radians and negate

        # Calculate derived constants
        self.CI = math.sqrt(self.GAM * self.AR * self.TEMP) # Check units of AR, TEMP
        print(f"K STAGE SPEED OF SOUND AT INLET  {self.CI:9.2f}")
        self.CP = self.AR * self.GAM / (self.GAM - 1.0)
        self.EXPON = 1.0 / (self.GAM - 1.0)

        # Precompute thickness interpolator
        # Ensure XZ, XR are sorted
        sort_xz = np.argsort(self.XZ)
        self.XZ = self.XZ[sort_xz]
        self.TN = self.TN[sort_xz, :]
        sort_xr = np.argsort(self.XR)
        self.XR = self.XR[sort_xr]
        self.TN = self.TN[:, sort_xr]

        self.tn_interpolator = RegularGridInterpolator(
            (self.XZ, self.XR), self.TN, method='linear', bounds_error=False, fill_value=0.0
        )

        print("Finished reading and initializing.")


    def _allocate_arrays(self):
        """Allocates NumPy arrays based on grid dimensions."""
        shape_mx_kmx = (self.MX, self.KMX)
        shape_mx = (self.MX,)
        shape_kmx = (self.KMX,)

        self.AL = np.zeros(shape_mx_kmx)
        self.BETA = np.zeros(shape_mx_kmx)
        self.CAL = np.zeros(shape_mx_kmx)
        self.CBETA = np.zeros(shape_mx_kmx)
        self.CURV = np.zeros(shape_mx_kmx)
        self.DN = np.zeros(shape_mx_kmx)
        self.PRS = np.zeros(shape_mx_kmx)
        self.R = np.zeros(shape_mx_kmx)
        self.Z = np.zeros(shape_mx_kmx)
        self.SM = np.zeros(shape_mx_kmx)
        self.SA = np.zeros(shape_mx_kmx)
        self.SB = np.zeros(shape_mx_kmx)
        self.SC = np.zeros(shape_mx_kmx)
        self.SD = np.zeros(shape_mx_kmx)
        self.SAL = np.zeros(shape_mx_kmx)
        self.SBETA = np.zeros(shape_mx_kmx)
        self.TT = np.zeros(shape_mx_kmx)
        self.WA = np.zeros(shape_mx_kmx)
        self.WTR = np.zeros(shape_mx_kmx)

        self.AB = np.zeros(shape_mx)
        self.AC = np.zeros(shape_mx)
        self.AD = np.zeros(shape_mx) # Also used for KMX size later
        self.DELBTA = np.zeros(shape_mx)
        self.DRDM = np.zeros(shape_mx)
        self.DTDR = np.zeros(shape_mx)
        self.DTDZ = np.zeros(shape_mx)
        self.DWMDM = np.zeros(shape_mx)
        self.DWTDM = np.zeros(shape_mx)
        self.WTFL = np.zeros(shape_kmx) # Used for integration result

        # Reallocate AD if needed for KMX size
        self.AD_kmx = np.zeros(shape_kmx)


    def _spline_fit_and_deriv(self, x, y):
        """
        Fits a cubic spline and calculates derivatives and curvature.
        Corresponds to parts of the main loop using SPLINE.

        Args:
            x (np.ndarray): Independent variable (e.g., transformed Z-R).
            y (np.ndarray): Dependent variable (e.g., transformed Z+R).

        Returns:
            tuple: (slope, curvature) arrays.
        """
        # Use natural boundary conditions, similar to the Fortran code's end conditions
        cs = CubicSpline(x, y, bc_type='natural')
        slope_raw = cs.derivative(nu=1)(x)
        deriv2 = cs.derivative(nu=2)(x)
        # Curvature = y'' / (1 + y'^2)^(3/2)
        curvature = deriv2 / (1.0 + slope_raw**2)**1.5
        return slope_raw, curvature

    def _spline_deriv_at_points(self, x_known, y_known, x_eval):
        """
        Calculates the derivative of a spline at specified points.
        Corresponds to SPLDER.

        Args:
            x_known (np.ndarray): X coordinates of known points.
            y_known (np.ndarray): Y coordinates of known points.
            x_eval (np.ndarray): X coordinates where derivative is needed.

        Returns:
            np.ndarray: Derivative values at x_eval.
        """
        # Ensure known points are sorted
        sort_idx = np.argsort(x_known)
        x_known_sorted = x_known[sort_idx]
        y_known_sorted = y_known[sort_idx]

        cs = CubicSpline(x_known_sorted, y_known_sorted, bc_type='natural')
        return cs.derivative(nu=1)(x_eval)

    def _spline_interpolate(self, x_known, y_known, x_eval):
        """
        Interpolates using a spline at specified points.
        Corresponds to SPLINT.

        Args:
            x_known (np.ndarray): X coordinates of known points.
            y_known (np.ndarray): Y coordinates of known points.
            x_eval (np.ndarray): X coordinates where interpolation is needed.

        Returns:
            np.ndarray: Interpolated values at x_eval.
        """
         # Ensure known points are sorted
        sort_idx = np.argsort(x_known)
        x_known_sorted = x_known[sort_idx]
        y_known_sorted = y_known[sort_idx]

        cs = CubicSpline(x_known_sorted, y_known_sorted, bc_type='natural')
        return cs(x_eval)

    def _cumulative_integral(self, x, y):
        """
        Calculates the cumulative integral using the trapezoidal rule.
        Corresponds to INTGRL (which used spline integration).
        Using trapezoidal rule for simplicity as SciPy spline doesn't
        directly give cumulative integral array easily.

        Args:
            x (np.ndarray): Independent variable.
            y (np.ndarray): Dependent variable.

        Returns:
            np.ndarray: Cumulative integral values.
        """
        if len(x) < 2:
            return np.zeros_like(x)
        # cumulative_trapezoid returns N-1 values, prepend 0
        integral = cumulative_trapezoid(y, x, initial=0)
        return integral

    def _linint(self, x1, y1):
        """
        Performs 2D linear interpolation on the blade thickness grid (TN).
        Corresponds to LININT.

        Args:
            x1 (float): Z coordinate for interpolation.
            y1 (float): R coordinate for interpolation.

        Returns:
            float: Interpolated thickness value.
        """
        # Uses the precomputed RegularGridInterpolator
        point = np.array([[x1, y1]])
        interpolated_value = self.tn_interpolator(point)[0]
        return interpolated_value

    def _contin(self, wa_current, wtfl_calc, ind_val, i_stream):
        """
        Adjusts WA(i, 0) based on calculated vs target weight flow.
        Corresponds to CONTIN subroutine logic. This is a simplified
        translation focusing on the core adjustment logic. The original
        had more complex state management (SPEED, WEIGHT arrays).

        Args:
            wa_current (float): Current WA[i, 0] value.
            wtfl_calc (float): Calculated weight flow for streamline i.
            ind_val (int): State indicator (not fully replicated).
            i_stream (int): Streamline index.

        Returns:
            tuple: (new_wa, new_ind) - Adjusted WA and next state indicator.
                   Returns None if convergence fails (original wrote error).
        """
        # Simplified logic: Proportional adjustment
        # The original Fortran logic is a custom iterative scheme.
        # A robust implementation might use scipy.optimize.root_scalar.
        # This is a basic proportional feedback:
        if abs(wtfl_calc - self.WT) <= self.WTOLER:
            return wa_current, 0 # Converged state for this streamline

        if wtfl_calc == 0: # Avoid division by zero
             # If flow is zero, maybe increase velocity slightly? Needs tuning.
             new_wa = wa_current * 1.1
        else:
             # Proportional adjustment factor
             # Clamp the adjustment factor to avoid extreme changes
             factor = max(0.5, min(2.0, self.WT / wtfl_calc))
             new_wa = wa_current * factor

        # The original IND logic managed a bracketing/search process.
        # We just return a state indicating iteration should continue.
        new_ind = 1 # Simplified: always indicate to continue iteration

        # Add a check to prevent WA from becoming excessively large or small/negative
        new_wa = max(1e-3, new_wa) # Prevent zero or negative velocity

        # Check for potential divergence (very basic)
        if new_wa > wa_current * 5 or new_wa < wa_current * 0.2:
             print(f"Warning: Large change in WA for streamline {i_stream}. "
                   f"Old={wa_current:.4f}, New={new_wa:.4f}, WTFL={wtfl_calc:.4f}")
             # Consider adding logic from original CONTIN for bracketing if needed

        # Original code had error exit (IND=6)
        # if some_divergence_condition:
        #    print(f"ERROR: Convergence failed for streamline {i_stream}, Max WT = {wtfl_calc}")
        #    return None, 6 # Signal failure

        return new_wa, new_ind


    def run_solver(self):
        """Performs the main iterative calculation loop."""

        initial_iter = self.ITER # Store initial iteration count setting

        while self.ITER >= 0: # Loop continues as long as ITER is non-negative
            if initial_iter == 0: # Only print headers if it was initially set to 0 (run till convergence)
                 print(f"\nIteration No. {self.itno:3d}")
                 print(f"{' ':6}{'AL':>9}{'RC':>9}{'SM':>9}{'BETA':>9}{'TT':>9}"
                       f"{'SA':>9}{'SB':>9}{'SC':>9}{'SD':>9}")

            self.errori = self.error # Store error from previous iteration
            self.error = 0.0        # Reset max error for current iteration

            # --- Calculate Geometric Parameters and Coefficients (Loops 160-230) ---
            for k in range(self.KMX): # Loop over streamlines
                # Transform coordinates for spline fitting
                # AB = (Z - R) / sqrt(2)
                # AC = (Z + R) / sqrt(2)
                # Using SM as the independent variable for splines along streamline
                # First calculate SM
                if k > 0: # SM depends on previous streamline's geometry (or initial guess)
                    # Recalculate SM based on current Z, R for this streamline k
                    self.SM[0, k] = 0.0
                    for i in range(1, self.MX):
                         ds = math.sqrt((self.Z[i, k] - self.Z[i-1, k])**2 +
                                        (self.R[i, k] - self.R[i-1, k])**2)
                         self.SM[i, k] = self.SM[i-1, k] + ds
                else: # k == 0 (hub streamline) - SM might be initialized differently or based on hub geometry
                      # Assuming SM is calculated based on Z, R of the hub itself if needed
                      # For now, keep SM[i, 0] as potentially zero or initialized elsewhere
                      pass # SM[i,0] might not be used directly in calcs if WA[i,0] is handled specially

                # Fit spline along the streamline (using SM as independent variable)
                # Need Z(SM) and R(SM) splines to get derivatives dZ/dSM, dR/dSM
                if self.SM[-1, k] > 1e-9: # Ensure SM has non-zero length
                    spline_z = CubicSpline(self.SM[:, k], self.Z[:, k], bc_type='natural')
                    spline_r = CubicSpline(self.SM[:, k], self.R[:, k], bc_type='natural')
                    dZdSM = spline_z.derivative(nu=1)(self.SM[:, k])
                    dRdSM = spline_r.derivative(nu=1)(self.SM[:, k])
                    d2ZdSM2 = spline_z.derivative(nu=2)(self.SM[:, k])
                    d2RdSM2 = spline_r.derivative(nu=2)(self.SM[:, k])

                    # Meridional angle AL = atan(dR/dZ) = atan( (dR/dSM) / (dZ/dSM) )
                    # Avoid division by zero if dZ/dSM is small
                    self.AL[:, k] = np.arctan2(dRdSM, dZdSM) # More robust than atan(dR/dZ)

                    # Curvature = (z'r'' - z''r') / (z'^2 + r'^2)^(3/2) where ' = d/dSM
                    denom = (dZdSM**2 + dRdSM**2)**1.5
                    # Avoid division by zero
                    denom[denom < 1e-9] = 1e-9
                    self.CURV[:, k] = (dZdSM * d2RdSM2 - d2ZdSM2 * dRdSM) / denom

                else: # Handle zero-length streamline (e.g., first point)
                    self.AL[:, k].fill(0.0) # Or estimate from geometry
                    self.CURV[:, k].fill(0.0) # Curvature is ill-defined

                self.CAL[:, k] = np.cos(self.AL[:, k])
                self.SAL[:, k] = np.sin(self.AL[:, k])

                # Calculate blade angle related terms (Loop 220)
                # Need dTdZ (derivative of thickness wrt Z along streamline)
                # Need dTdR (derivative of thickness wrt R along streamline)

                # Use spline for THTA vs XT (blade angle vs thickness coord?) - Fortran uses Z(I,K)
                # Assuming XT are coordinates along Z where THTA is defined
                # Check if XT corresponds to Z or SM
                # Fortran uses CALL SPLDER(XT(1), THTA(1), MTHTA, Z(1,K), MX, DTDZ(1))
                 # This implies THTA is a function of Z, defined by points (XT, THTA)
                # Let's assume XT are Z-coordinates
                # self.DTDZ[:,k] = self._spline_deriv_at_points(self.XT, self.THTA, self.Z[:, k]) # Incorrect indexing
                self.DTDZ = self._spline_deriv_at_points(self.XT, self.THTA, self.Z[:, k]) # Corrected: DTDZ is 1D (size MX)

                # Calculate dTdR based on formula involving RB
                rinlet = (self.RS[0] + self.RH[0]) / 2.0 # Inlet radius (approx)
                # Calculate dTdR based on formula involving RB
                rinlet = (self.RS[0] + self.RH[0]) / 2.0 # Inlet radius (approx)
                # CEF = SIN(BETIN)/COS(BETIN)/RINLET/(RINLET-RB)**2 - Check BETIN sign usage
                cef_term = math.tan(self.BETIN) / rinlet / (rinlet - self.RB)**2 if abs(rinlet - self.RB) > 1e-9 else 0.0

                for i in range(self.MX):
                    # Interpolate thickness T at current (Z, R)
                    # T = self._linint(self.Z[i, k], self.R[i, k]) # Original Fortran call
                    # The Fortran code seems to use T for TT calculation, not directly for dTdR/dTdZ here.
                    # Let's follow the Fortran logic for DTDR and DTDZ calculation first.

                    if self.R[i, k] <= self.RB:
                        self.DTDR[i] = 0.0
                    else:
                        # DTDR(I) = CEF*(R(I,K)-RB)**2 - This seems like d(Theta)/dR, not d(Thickness)/dR
                        # Let's assume DTDR is dTheta/dR based on context
                        self.DTDR[i] = cef_term * (self.R[i, k] - self.RB)**2

                    # Interpolate blade thickness T using LININT
                    t_thick = self._linint(self.Z[i, k], self.R[i, k])

                    # Calculate TT (Tangential thickness) and BETA (Relative flow angle)
                    tq = self.R[i, k] * self.DTDR[i] # R * dTheta/dR
                    tp = self.R[i, k] * self.DTDZ[i] # R * dTheta/dZ (Theta from spline)

                    # TT(I,K) = T*SQRT(1.0+TP*TP) - Check if T is thickness or angle Theta
                    # Assuming T is thickness from LININT
                    self.TT[i, k] = t_thick * math.sqrt(1.0 + tp**2) # Tangential thickness

                    # BETA(I,K) = ATAN(TP*CAL(I,K)+TQ*SAL(I,K))
                    self.BETA[i, k] = math.atan(tp * self.CAL[i, k] + tq * self.SAL[i, k])
                    self.SBETA[i, k] = math.sin(self.BETA[i, k])
                    self.CBETA[i, k] = math.cos(self.BETA[i, k])

                    # Calculate coefficients SA, SC (depend on geometry and flow angle)
                    self.SA[i, k] = (self.CBETA[i, k]**2 * self.CAL[i, k] * self.CURV[i, k]
                                     - self.SBETA[i, k]**2 / self.R[i, k] if self.R[i, k] > 1e-9 else 0.0
                                     + self.SAL[i, k] * self.CBETA[i, k] * self.SBETA[i, k] * self.DTDR[i]) # Using DTDR = dTheta/dR

                    self.SC[i, k] = (-self.SAL[i, k] * self.CBETA[i, k]**2 * self.CURV[i, k]
                                     + self.SAL[i, k] * self.CBETA[i, k] * self.SBETA[i, k] * self.DTDZ[i]) # Using DTDZ = dTheta/dZ

                # Calculate velocity derivatives using splines (DWMDM, DWTDM)
                # AB = WA * CBETA (Meridional velocity component?)
                # AC = WA * SBETA (Tangential velocity component relative to meridional?)
                self.AB[:] = self.WA[:, k] * self.CBETA[:, k]
                self.AC[:] = self.WA[:, k] * self.SBETA[:, k]

                if self.SM[-1, k] > 1e-9: # Check if streamline has length
                    spline_AB = CubicSpline(self.SM[:, k], self.AB, bc_type='natural')
                    spline_AC = CubicSpline(self.SM[:, k], self.AC, bc_type='natural')
                    self.DWMDM[:] = spline_AB.derivative(nu=1)(self.SM[:, k]) # d(W_meridional)/d(SM)
                    self.DWTDM[:] = spline_AC.derivative(nu=1)(self.SM[:, k]) # d(W_tangential_rel)/d(SM)
                else:
                    self.DWMDM[:].fill(0.0)
                    self.DWTDM[:].fill(0.0)

                # Calculate coefficients SB, SD (depend on velocity derivatives)
                for i in range(self.MX):
                    # SB(I,K) = SAL(I,K)*CBETA(I,K)*DWMDM(I)-2.0*W*SBETA(I,K)
                    #    1    +DTDR(I)*R(I,K)*CBETA(I,K)*(DWTDM(I)+2.0*W*SAL(I,K))
                    self.SB[i, k] = (self.SAL[i, k] * self.CBETA[i, k] * self.DWMDM[i]
                                     - 2.0 * self.W * self.SBETA[i, k]
                                     + self.DTDR[i] * self.R[i, k] * self.CBETA[i, k] *
                                       (self.DWTDM[i] + 2.0 * self.W * self.SAL[i, k]))

                    # SD(I,K) = CAL(I,K)*CBETA(I,K)*DWMDM(I)+DTDZ(I)*
                    #    1    R(I,K)*CBETA(I,K)*(DWTDM(I)+2.0*W*SAL(I,K))
                    self.SD[i, k] = (self.CAL[i, k] * self.CBETA[i, k] * self.DWMDM[i]
                                     + self.DTDZ[i] * self.R[i, k] * self.CBETA[i, k] *
                                       (self.DWTDM[i] + 2.0 * self.W * self.SAL[i, k]))

                    # --- Print intermediate results (if initial ITER <= 0) ---
                    if initial_iter <= 0 and k % self.NPRT == 0:
                        if i == 0: # Print streamline header only once
                             print(f"  Streamline{k+1:3d}")
                        a_deg = self.AL[i, k] * self.RAD_TO_DEG
                        b_sm_in = self.SM[i, k] / self.IN_TO_FT
                        e_tt_in = self.TT[i, k] / self.IN_TO_FT
                        g_beta_deg = self.BETA[i, k] * self.RAD_TO_DEG
                        print(f"{a_deg:14.6f}{self.CURV[i, k]:14.6f}{b_sm_in:14.6f}"
                              f"{g_beta_deg:14.6f}{e_tt_in:14.6f}"
                              f"{self.SA[i, k]:14.6f}{self.SB[i, k]:14.6f}"
                              f"{self.SC[i, k]:14.6f}{self.SD[i, k]:14.6f}")

            # --- Calculate Blade Surface Velocities (WTR) (Loop 250) ---
            # This is done *after* convergence in the Fortran code (ITER == 0)
            # Here we calculate it every iteration if initial_iter <= 0
            if initial_iter <= 0:
                for k in range(self.KMX):
                    if self.SM[-1, k] > 1e-9: # Check streamline length
                        # DELBTA = d(TT)/d(SM)
                        spline_TT = CubicSpline(self.SM[:, k], self.TT[:, k], bc_type='natural')
                        self.DELBTA[:] = spline_TT.derivative(nu=1)(self.SM[:, k])

                        # AB = (R*W + WA*SBETA) * (Circumference / N - TT)
                        # Need DRDM = d(AB)/d(SM)
                        circum_factor = 2.0 * self.PI / self.XN
                        self.AB[:] = ((self.R[:, k] * self.W + self.WA[:, k] * self.SBETA[:, k]) *
                                      (circum_factor * self.R[:, k] - self.TT[:, k]))
                        spline_AB_surf = CubicSpline(self.SM[:, k], self.AB, bc_type='natural')
                        self.DRDM[:] = spline_AB_surf.derivative(nu=1)(self.SM[:, k]) # Store in DRDM

                        # Handle splitter factor SFACT > 1.0
                        ad_temp = np.zeros_like(self.DRDM) # Temporary array like AD in Fortran
                        if self.SFACT > 1.0:
                            circum_factor_s = 2.0 * self.PI / (self.SFACT * self.XN)
                            self.AB[:] = ((self.R[:, k] * self.W + self.WA[:, k] * self.SBETA[:, k]) *
                                          (circum_factor_s * self.R[:, k] - self.TT[:, k]))
                            spline_AB_surf_s = CubicSpline(self.SM[:, k], self.AB, bc_type='natural')
                            ad_temp[:] = spline_AB_surf_s.derivative(nu=1)(self.SM[:, k]) # Store in temp AD

                        for i in range(self.MX):
                            betad = self.BETA[i, k] - self.DELBTA[i] / 2.0 # Pressure surface angle
                            betat = betad + self.DELBTA[i] # Suction surface angle
                            cosbd = math.cos(betad)
                            cosbt = math.cos(betat)

                            # Use AD if Z < ZSPLIT and SFACT > 1.0
                            drdm_eff = ad_temp[i] if self.SFACT > 1.0 and self.Z[i, k] < self.ZSPLIT else self.DRDM[i]

                            # Avoid division by zero for cosbd + cosbt
                            denom_wtr = cosbd + cosbt
                            if abs(denom_wtr) < 1e-9:
                                self.WTR[i, k] = 0.0 # Or handle appropriately
                            else:
                                # WTR(I,K) = COSBD*COSBT/(COSBD+COSBT)*(2.0*WA(I,K)/COSBD+R(I,K)
                                #    1    *W*(BETAD-BETAT)/CBETA(I,K)**2+DRDM(I))
                                term1 = 2.0 * self.WA[i, k] / cosbd if abs(cosbd) > 1e-9 else 0.0
                                term2 = (self.R[i, k] * self.W * (betad - betat) / self.CBETA[i, k]**2
                                         if abs(self.CBETA[i, k]) > 1e-9 else 0.0)
                                self.WTR[i, k] = (cosbd * cosbt / denom_wtr *
                                                  (term1 + term2 + drdm_eff))
                    else:
                        self.WTR[:, k].fill(0.0) # No velocity if streamline has no length


            # --- Calculate Weight Flow and Update Streamlines (Loops 260-370) ---
            for i in range(self.MX): # Loop over quasi-orthogonals
                ind = 1 # Initial state for CONTIN logic
                wa_i0_current = self.WA[i, 0] # Store initial WA at hub for this orthogonal

                # Inner loop to converge WA[i,0] to match WT using CONTIN logic
                contin_attempts = 0
                max_contin_attempts = 20 # Limit attempts to prevent infinite loops

                while ind != 0 and contin_attempts < max_contin_attempts:
                    contin_attempts += 1
                    # Update WA along the quasi-orthogonal (k=1 to KMX-1) based on WA[i,0] and coeffs
                    # WA(I,K) = (WAS+WASS)/2.0 (Predictor-corrector?)
                    for k in range(1, self.KMX):
                        j = k - 1
                        hr = self.R[i, k] - self.R[i, j]
                        hz = self.Z[i, k] - self.Z[i, j]

                        # Predictor step (WAS)
                        was = self.WA[i, j] * (1.0 + self.SA[i, j] * hr + self.SC[i, j] * hz)

                        # Corrector step uses predicted WAS to estimate influence at k
                        # WASS = WA(I,J)+WAS*(SA(I,K)*HR+SC(I,K)*HZ) + SB(I,K)*HR+SD(I,K)*HZ
                        # Note: Fortran used WA(I,J) in WASS, seems like typo, should be WAS?
                        # Let's assume it used WA(I,J) as written.
                        wass = (self.WA[i, j] +
                                was * (self.SA[i, k] * hr + self.SC[i, k] * hz) +
                                self.SB[i, k] * hr + self.SD[i, k] * hz)

                        # Average predictor and corrector
                        self.WA[i, k] = (was + wass) / 2.0

                        # Ensure WA doesn't become negative (physical constraint)
                        if self.WA[i, k] < 0:
                             # print(f"Warning: Negative WA[{i},{k}] calculated. Setting to small positive.")
                             self.WA[i, k] = 1e-3


                    # Calculate density, pressure, and integrand for weight flow (Loop 340)
                    self.AD_kmx.fill(0.0) # Reset integrand array (size KMX)
                    for k in range(self.KMX):
                        # Thermodynamics: Calculate Temperature ratio TIP
                        # TIP = 1.0-(WA(I,K)**2+2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
                        # Ensure denominator is not zero
                        denom_tip = 2.0 * self.CP * self.TEMP
                        if abs(denom_tip) < 1e-9:
                            tip = 1.0 # Avoid division by zero, maybe handle error?
                        else:
                            tip = 1.0 - (self.WA[i, k]**2 + 2.0 * self.W * self.ALM - (self.W * self.R[i, k])**2) / denom_tip

                        # Check for non-physical TIP (negative absolute temperature)
                        if tip < 0.0:
                            # Fortran jumps back to 280 (WA(I,1) = 0.5*WA(I,1))
                            # print(f"Warning: Negative TIP encountered at i={i}, k={k}. Halving WA[{i},0].")
                            wa_i0_current *= 0.5
                            self.WA[i, 0] = wa_i0_current # Adjust hub velocity
                            ind = 1 # Force re-calculation of WA along orthogonal
                            break # Exit k loop, restart inner while loop with new WA[i,0]

                        # Calculate TPPIP (related to total conditions?)
                        # TPPIP = 1.0-(2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
                        tppip = 1.0 - (2.0 * self.W * self.ALM - (self.W * self.R[i, k])**2) / denom_tip if abs(denom_tip) > 1e-9 else 1.0

                        # Calculate density with loss term
                        # DENSTY = TIP**EXPON*RHO-(TIP/TPPIP)**EXPON*
                        #    1    PLOSS/AR/TPPIP/TEMP*32.17*SM(I,K)/SM(MX,K)
                        # Need SM(MX,K) - meridional length of the k-th streamline
                        sm_mx_k = self.SM[self.MX-1, k] if self.MX > 0 else 0.0
                        loss_term = 0.0
                        if abs(self.AR * tppip * self.TEMP) > 1e-9 and sm_mx_k > 1e-9 and tppip > 1e-9:
                            loss_term = ((tip / tppip)**self.EXPON * self.PLOSS /
                                         (self.AR * tppip * self.TEMP) * 32.17 * # Assuming 32.17 is g
                                         self.SM[i, k] / sm_mx_k)

                        density = tip**self.EXPON * self.RHO - loss_term
                        if density < 0:
                            # print(f"Warning: Negative density calculated at i={i}, k={k}. Setting to small positive.")
                            density = 1e-4

                        # Calculate pressure (PRS)
                        # PRS(I,K) = DENSTY*AR*TIP*TEMP/32.17/144.0
                        self.PRS[i, k] = density * self.AR * tip * self.TEMP / 32.17 / self.PSI_TO_PSF # Convert psf to psi

                        # Calculate integrand for weight flow (AD(K))
                        # PSI = angle of the quasi-orthogonal
                        if abs(self.ZS[i] - self.ZH[i]) > 1e-9:
                            psi = math.atan((self.RS[i] - self.RH[i]) / (self.ZS[i] - self.ZH[i])) - self.PI / 2.0
                        elif abs(self.RS[i] - self.RH[i]) > 1e-9:
                             # Handle vertical orthogonal line
                             psi = 0.0 if self.RS[i] > self.RH[i] else self.PI # Check direction
                        else:
                             psi = 0.0 # Degenerate case

                        # WTHRU = WA(I,K)*CBETA(I,K)*COS(PSI-AL(I,K)) (Velocity normal to orthogonal)
                        wthru = self.WA[i, k] * self.CBETA[i, k] * math.cos(psi - self.AL[i, k])

                        # Effective number of blades (consider splitter)
                        a_blades = self.SFACT * self.XN if self.Z[i, k] < self.ZSPLIT else self.XN

                        # C = passage width = 2*pi*R - N*TT
                        c_width = 2.0 * self.PI * self.R[i, k] - a_blades * self.TT[i, k]
                        if c_width < 0:
                            # print(f"Warning: Negative passage width C at i={i}, k={k}. Clamping.")
                            c_width = 1e-4 # Avoid negative flow

                        # Integrand AD(K) = density * W_normal * Area_element
                        # Area element is C * d(DN)
                        self.AD_kmx[k] = density * wthru * c_width

                    else: # This else corresponds to the 'for k in range(self.KMX)' loop
                        # If the inner loop finished without break (no negative TIP)
                        # Integrate AD vs DN to get weight flow
                        # AC(K) = DN(I,K) in Fortran
                        self.WTFL = self._cumulative_integral(self.DN[i, :], self.AD_kmx)
                        wtfl_calc = self.WTFL[-1] # Total calculated weight flow for this orthogonal

                        # Check convergence and adjust WA[i, 0] using CONTIN logic
                        if abs(self.WT - wtfl_calc) <= self.WTOLER:
                            ind = 0 # Converged for this orthogonal
                        else:
                            # Call the CONTIN equivalent
                            wa_i0_current, ind = self._contin(wa_i0_current, wtfl_calc, ind, i)
                            if ind == 6: # Check if CONTIN signalled failure
                                print(f"Error: CONTIN failed for streamline {i}")
                                self.ITER = -1 # Force exit
                                break # Exit inner while loop
                            self.WA[i, 0] = wa_i0_current # Update hub velocity for next attempt

                        if ind != 0 and contin_attempts >= max_contin_attempts:
                             print(f"Warning: Max CONTIN attempts reached for orthogonal {i}. "
                                   f"Calculated WTFL={wtfl_calc:.4f}, Target WT={self.WT:.4f}")
                             ind = 0 # Stop trying for this orthogonal

                # --- Update Streamline Positions (DN) ---
                if self.ITER >= 0: # Only update if not exiting due to error
                    # Interpolate DN based on target weight flow distribution BA
                    # Fortran uses SPLINT(WTFL, AC, KMX, BA, KMX, AB) where AC=DN, BA=Target WT dist.
                    # Result AB is the DN corresponding to the target BA distribution.
                    # Ensure WTFL is monotonically increasing for interpolation
                    if np.all(np.diff(self.WTFL) >= 0):
                         # Use linear interpolation as spline might be unstable if WTFL is noisy
                         # interp_func = CubicSpline(self.WTFL, self.DN[i, :], bc_type='natural')
                         # dn_new = interp_func(self.BA)
                         dn_new = np.interp(self.BA, self.WTFL, self.DN[i, :])
                    else:
                         # Handle non-monotonic WTFL - maybe keep old DN?
                         print(f"Warning: Non-monotonic WTFL encountered for orthogonal {i}. Skipping DN update.")
                         dn_new = self.DN[i, :].copy() # Keep old values


                    # Calculate change and update DN with correction factor
                    delta = np.abs(dn_new - self.DN[i, :])
                    self.DN[i, :] = (1.0 - self.CORFAC) * self.DN[i, :] + self.CORFAC * dn_new

                    # Update maximum error
                    max_delta_i = np.max(delta) if delta.size > 0 else 0.0
                    if max_delta_i > self.error:
                        self.error = max_delta_i

            # End of loop for quasi-orthogonals (i)

            # --- Update Streamline Coordinates Z, R (Loop 380) ---
            # Update interior streamlines based on new DN distribution
            for k in range(1, self.KMXM1): # Loop from k=1 to KMX-2
                for i in range(self.MX):
                    # Linear interpolation between hub (k=0) and shroud (k=KMX-1) geometry
                    # based on the ratio of normal distances DN
                    dn_ratio = self.DN[i, k] / self.DN[i, self.KMXM1] if abs(self.DN[i, self.KMXM1]) > 1e-9 else 0.0
                    # Clamp ratio to avoid extrapolation issues if DN is weird
                    dn_ratio = max(0.0, min(1.0, dn_ratio))

                    self.Z[i, k] = dn_ratio * (self.ZS[i] - self.ZH[i]) + self.ZH[i]
                    self.R[i, k] = dn_ratio * (self.RS[i] - self.RH[i]) + self.RH[i]

            # --- Check Convergence and Decrement Iteration Counter ---
            if self.error <= self.TOLER:
                print(f"Converged after {self.itno} iterations.")
                self.ITER = -1 # Stop iteration (like Fortran setting ITER < 0)
            elif self.error >= self.errori and self.itno > 1: # Check if error increased (divergence)
                print(f"Warning: Error increased from {self.errori:.6f} to {self.error:.6f}. Stopping.")
                self.ITER = -1 # Stop iteration
            else:
                # Normal iteration decrement (only if initial ITER > 0)
                if initial_iter > 0:
                    self.ITER -= 1

                # If initial_iter was 0, we run until convergence or divergence,
                # so don't decrement ITER based on its value.
                # The loop condition `while self.ITER >= 0:` handles termination.

            # Increment iteration number for next loop or final output
            self.itno += 1

            # Safety break if initial_iter was 0 and it runs too long
            if initial_iter == 0 and self.itno > 500: # Arbitrary limit
                 print("Warning: Maximum iterations (500) reached without convergence.")
                 self.ITER = -1 # Force stop


        # --- End of Iteration Loop ---
        print("\n--- Final Results ---")
        self.write_output()

        # Final error print
        print(f"\nIteration No. {self.itno-1:3d}    Max. Streamline Change = {self.error / self.IN_TO_FT:10.6f} (inches)")


    def write_output(self):
        """Prints the final calculated results."""
        print(f"\n{' ':9}{'DN':>19}{'Z':>19}{'R':>19}"
              f"{'WA':>19}{'PRS':>19}{'WTR':>18}{'RC':>18}")

        # Recalculate final geometry parameters (AL, CURV) if needed for output
        # (They might be slightly different after the last DN update)
        # Optional: Add recalculation of AL, CURV here if exact final values are needed

        for k in range(0, self.KMX, self.NPRT): # Print every NPRT streamline
            print(f"\n  Streamline{k+1:3d}")
            for i in range(self.MX):
                # Convert back to inches and psi for printing if desired
                dn_pr = self.DN[i, k] / self.IN_TO_FT
                z_pr = self.Z[i, k] / self.IN_TO_FT
                r_pr = self.R[i, k] / self.IN_TO_FT
                prs_pr = self.PRS[i, k] # Already in PSI
                # WTR and WA are velocities (ft/s)
                # RC is curvature (1/ft)

                # Need to recalculate curvature for the final geometry if not done after loop
                # Placeholder for curvature - using last calculated value
                curv_pr = self.CURV[i, k] # Might be slightly off

                print(f"{dn_pr:19.6f}{z_pr:19.6f}{r_pr:19.6f}"
                      f"{self.WA[i, k]:19.6f}{prs_pr:19.6f}"
                      f"{self.WTR[i, k]:18.6f}{curv_pr:18.6f}")

    def dump_bcd(self):
        """Dumps final arrays in a format readable by BCREAD (placeholder)."""
        if self.BCDP != 1:
            return
        print("\n--- Dumping BCD Data ---")
        # This needs to format the output exactly as BCREAD expects.
        # Assuming simple space-separated floats per line.
        print("DN array:")
        for i in range(self.MX):
            print(*(f"{x:15.8E}" for x in self.DN[i, :])) # Example format
        print("WA array:")
        for i in range(self.MX):
            print(*(f"{x:15.8E}" for x in self.WA[i, :]))
        print("Z array:")
        for i in range(self.MX):
            print(*(f"{x:15.8E}" for x in self.Z[i, :]))
        print("R array:")
        for i in range(self.MX):
            print(*(f"{x:15.8E}" for x in self.R[i, :]))
        print("--- End BCD Dump ---")

# --- Main execution ---
if __name__ == "__main__":
    solver = TurbomachinerySolver()
    try:
        # Loop to handle multiple runs like the original GOTO 10
        while True:
            try:
                solver.read_input()
                solver.run_solver()
                solver.dump_bcd()
                # Check if more input is available (e.g., for another run)
                # This part depends on how input is provided.
                # If reading from stdin, attempting another read might work,
                # or expect a specific signal/empty line to stop.
                # For simplicity, we break after one run here.
                # To enable multiple runs, you might need to check for EOF
                # or specific input termination.
                print("\nRun complete. Add logic to check for more input data for next run.")
                break # Exit after the first successful run

            except EOFError:
                print("\nEnd of input reached.")
                break # Exit loop if no more input
            except Exception as e:
                print(f"\nAn error occurred during processing: {e}")
                import traceback
                traceback.print_exc()
                # Decide whether to stop or try reading next dataset if possible
                break # Stop on error for now

    except KeyboardInterrupt:
        print("\nExecution interrupted by user.")
