import numpy as np
from scipy.interpolate import CubicSpline # Removed RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.special import comb # For Bezier binomial coefficients
import math
import matplotlib.pyplot as plt

# --- Bezier Curve Evaluation ---
def evaluate_bezier(control_points, num_points):
    """
    Evaluates a Bezier curve defined by control_points at num_points intervals.

    Args:
        control_points (np.ndarray): Array of control points, shape (n+1, dim).
        num_points (int): Number of points to evaluate on the curve.

    Returns:
        np.ndarray: Array of evaluated points, shape (num_points, dim).
    """
    n = control_points.shape[0] - 1 # Degree of the curve
    dim = control_points.shape[1]   # Dimension of the points (e.g., 2 for Z,R)
    t = np.linspace(0, 1, num_points) # Parameter t from 0 to 1

    curve_points = np.zeros((num_points, dim))
    for i in range(n + 1):
        # Bernstein polynomial basis function
        bernstein_poly = comb(n, i) * (1 - t)**(n - i) * t**i
        # Add contribution of control point P_i
        curve_points += np.outer(bernstein_poly, control_points[i])

    return curve_points

class TurbomachinerySolver:
    """
    Solves turbomachinery flow using a streamline curvature method.
    Uses Bezier curves for geometry and blade definition.
    """
    def __init__(self, mx_points=200): # Default number of QO lines
        """Initializes solver parameters and data arrays."""
        # Constants (Set some defaults, others in setup)
        self.RHO = 0.0
        self.GAM = 1.4
        self.AR = 1715.0 # Example: ft*lbf/(slug*R) for air
        self.TEMP = 592.0 # Example: Rankine
        self.W = 5390.0 * (2 * math.pi / 60.0) # Example: RPM to rad/s
        self.WT = 0.9840 # Example: Target weight flow (slug/s ?)
        self.ALM = 155.3 # Example: Inlet swirl parameter (r*Vt)
        self.TOLER = 0.001 # Example: Convergence tolerance (feet)
        self.PLOSS = 2.5 * 144.0 # Example: Pressure loss (psi to psf)
        self.WTOLER = 0.001 # Example: Convergence tolerance for weight flow
        self.NPRT = 2   # Print frequency for streamlines
        self.ITER = 0   # Run until convergence flag
        self.SFACT = 1.0 # Splitter factor
        self.ZSPLIT = -1.0 # Z-coordinate for splitter effect (inactive if < min Z)
        self.BETIN = -35.0 # Inlet blade angle (degrees)
        self.RB = 1.75 # Base radius for blade angle calculation (inches)
        self.CORFAC = 0.1 # Correction factor for streamline update
        self.XN = 13.0   # Number of blades
        self.KMX = 11    # Number of streamlines (hardcoded example)

        # Grid dimensions
        self.MX = mx_points # Number of points along meridional direction (QOs)
        # self.KMX = 0    # Set based on input/default
        # self.MR = 0     # No longer needed for TN grid
        # self.MZ = 0     # No longer needed for TN grid

        # Input geometry and blade data (will be generated from Bezier)
        self.ZS = None  # Shroud Z coordinates (MX)
        self.ZH = None  # Hub Z coordinates (MX)
        self.RS = None  # Shroud R coordinates (MX)
        self.RH = None  # Hub R coordinates (MX)
        # self.THTA = None # Blade angle data (MX) - Generated
        # self.XT = None  # Blade Z coordinate data (MX) - Generated
        self.blade_angle_spline = None # Spline for angle(Z)
        self.blade_thickness_spline = None # Spline for thickness(Z)
        # self.TN = None  # Blade thickness grid (MZ, MR) - Removed
        # self.XZ = None  # Axial coordinates for thickness grid (MZ) - Removed
        # self.XR = None  # Radial coordinates for thickness grid (MR) - Removed

        # Calculated variables (Arrays typically MX x KMX)
        self.AL = None
        self.BETA = None
        self.CAL = None
        self.CBETA = None
        self.CURV = None
        self.DN = None
        self.PRS = None
        self.R = None
        self.Z = None
        self.SM = None
        self.SA = None
        self.SB = None
        self.SC = None
        self.SD = None
        self.SAL = None
        self.SBETA = None
        # self.TN_interp = None # Removed
        self.TT = None
        self.WA = None
        self.WTR = None

        # Intermediate calculation arrays (typically size MX or KMX)
        self.AB = None
        self.AC = None
        self.AD = None
        self.BA = None
        self.DELBTA = None
        self.DRDM = None
        self.DTDR = None # dTheta/dR - Kept for radial variation effect
        self.DTDZ = None # dTheta/dZ - Calculated from angle spline
        self.DWMDM = None
        self.DWTDM = None
        self.WTFL = None

        # Control flags
        # self.TYPE = 0   # No longer needed
        self.BCDP = 0   # BCD dump flag (can keep if needed)
        # self.SRW = 0    # Removed

        # Calculated constants
        self.CI = 0.0
        self.CP = 0.0
        self.EXPON = 0.0
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

    def setup_from_bezier(self):
        """Sets up geometry, blade data, and parameters from hardcoded Bezier curves."""
        print("Setting up solver from Bezier curves...")
        self.runo += 1
        self.KMXM1 = self.KMX - 1

        # --- Hardcode Bezier Control Points (Z, R in inches) ---
        # Example: 5 control points for a simple converging-diverging nozzle shape
        # Hub: (Z, R) inches
        hub_cp = np.array([
            [0.0, 0.8],   # Inlet
            [0.3, 0.7],   # Converging
            [0.6, 0.65],  # Throat approx
            [0.9, 0.75],  # Diverging
            [1.2, 0.8]    # Outlet
        ]) * np.array([1.0, 1.0]) # Scale if needed

        # Shroud: (Z, R) inches
        shroud_cp = np.array([
            [0.0, 2.25],  # Inlet
            [0.3, 2.1],   # Converging
            [0.6, 2.0],   # Throat approx
            [0.9, 2.15],  # Diverging
            [1.2, 2.25]   # Outlet
        ]) * np.array([1.0, 1.0]) # Scale if needed

        # Blade Angle: (Z, Angle in degrees) - Angle relative to axial
        # Example: Starts negative (angled wrt axial), moves towards zero
        angle_cp = np.array([
            [0.0, -35.0], # Inlet Z, Inlet Angle (deg)
            [0.3, -25.0],
            [0.6, -10.0],
            [0.9, -5.0],
            [1.2, 0.0]   # Outlet Z, Outlet Angle (deg)
        ])

        # Blade Thickness: (Z, Thickness in inches)
        # Example: Simple thickness profile
        thickness_cp = np.array([
            [0.0, 0.08], # Inlet Z, Inlet Thickness (in)
            [0.3, 0.12],
            [0.6, 0.10], # Max thickness approx
            [0.9, 0.07],
            [1.2, 0.05]  # Outlet Z, Outlet Thickness (in)
        ])

        # --- Generate Geometry Points ---
        hub_points = evaluate_bezier(hub_cp, self.MX)
        shroud_points = evaluate_bezier(shroud_cp, self.MX)
        angle_points = evaluate_bezier(angle_cp, self.MX)
        thickness_points = evaluate_bezier(thickness_cp, self.MX)

        self.ZH = hub_points[:, 0]
        self.RH = hub_points[:, 1]
        self.ZS = shroud_points[:, 0]
        self.RS = shroud_points[:, 1]

        # Store generated angle/thickness data (Z coordinates might not be perfectly monotonic from Bezier)
        # Use the Z-coordinates from the angle curve evaluation as the basis for angle/thickness splines
        z_coords_blade = angle_points[:, 0]
        blade_angles_deg = angle_points[:, 1]
        # Interpolate thickness onto the same Z coordinates
        blade_thickness = np.interp(z_coords_blade, thickness_points[:, 0], thickness_points[:, 1])

        # --- Unit Conversions ---
        self.ZH *= self.IN_TO_FT
        self.RH *= self.IN_TO_FT
        self.ZS *= self.IN_TO_FT
        self.RS *= self.IN_TO_FT
        z_coords_blade *= self.IN_TO_FT
        blade_thickness *= self.IN_TO_FT
        blade_angles_rad = blade_angles_deg * self.DEG_TO_RAD # Keep angles in radians internally
        self.RB *= self.IN_TO_FT
        self.ZSPLIT *= self.IN_TO_FT
        self.TOLER *= self.IN_TO_FT
        self.BETIN *= self.DEG_TO_RAD # Convert input BETIN (used for DTDR)

        # --- Create Splines for Blade Angle and Thickness ---
        # Ensure Z coordinates are sorted for spline creation
        sort_idx = np.argsort(z_coords_blade)
        z_coords_blade_sorted = z_coords_blade[sort_idx]
        blade_angles_rad_sorted = blade_angles_rad[sort_idx]
        blade_thickness_sorted = blade_thickness[sort_idx]

        self.blade_angle_spline = CubicSpline(z_coords_blade_sorted, blade_angles_rad_sorted, bc_type='natural')
        self.blade_thickness_spline = CubicSpline(z_coords_blade_sorted, blade_thickness_sorted, bc_type='natural', extrapolate=False) # No extrapolation for thickness

        # --- Allocate Arrays ---
        self._allocate_arrays() # Uses self.MX, self.KMX

        # --- Initialize Flow Field ---
        # Similar to TYPE=0 logic, using generated geometry
        inlet_area_corrected = self.PI * (self.RS[0] + self.RH[0]) * (self.ZS[0] - self.ZH[0]) if abs(self.ZS[0] - self.ZH[0]) > 1e-9 else self.PI * (self.RS[0]**2 - self.RH[0]**2)
        wa_inlet = self.WT / (self.RHO * inlet_area_corrected) if self.RHO > 1e-9 and inlet_area_corrected > 1e-9 else 100.0 # Default guess if density/area is zero

        for i in range(self.MX):
            dn_kmx = math.sqrt((self.ZS[i] - self.ZH[i])**2 + (self.RS[i] - self.RH[i])**2)
            for k in range(self.KMX):
                frac = k / self.KMXM1 if self.KMXM1 > 0 else 0.5
                self.DN[i, k] = frac * dn_kmx
                self.WA[i, k] = wa_inlet # Initial guess
                self.Z[i, k] = frac * (self.ZS[i] - self.ZH[i]) + self.ZH[i]
                self.R[i, k] = frac * (self.RS[i] - self.RH[i]) + self.RH[i]
        self.DN[:, 0] = 0.0 # Hub streamline DN is 0

        # --- Initialize Other Parameters ---
        self.SM.fill(0.0)
        self.BA = np.linspace(0.0, self.WT, self.KMX) # Weight flow distribution guess

        # --- Calculate Derived Constants ---
        self.CI = math.sqrt(self.GAM * self.AR * self.TEMP) if self.TEMP > 0 else 0.0
        self.CP = self.AR * self.GAM / (self.GAM - 1.0) if self.GAM > 1.0 else 0.0
        self.EXPON = 1.0 / (self.GAM - 1.0) if self.GAM > 1.0 else 0.0

        print(f"Setup complete. MX={self.MX}, KMX={self.KMX}")
        print(f"Hub Z range: {self.ZH[0]/self.IN_TO_FT:.3f} to {self.ZH[-1]/self.IN_TO_FT:.3f} in")
        print(f"Shroud Z range: {self.ZS[0]/self.IN_TO_FT:.3f} to {self.ZS[-1]/self.IN_TO_FT:.3f} in")
        print(f"Hub R range: {self.RH[0]/self.IN_TO_FT:.3f} to {self.RH[-1]/self.IN_TO_FT:.3f} in")
        print(f"Shroud R range: {self.RS[0]/self.IN_TO_FT:.3f} to {self.RS[-1]/self.IN_TO_FT:.3f} in")
        print(f"Blade Angle Z range: {z_coords_blade_sorted[0]/self.IN_TO_FT:.3f} to {z_coords_blade_sorted[-1]/self.IN_TO_FT:.3f} in")
        print(f"Blade Angle range: {blade_angles_rad_sorted[0]*self.RAD_TO_DEG:.2f} to {blade_angles_rad_sorted[-1]*self.RAD_TO_DEG:.2f} deg")
        print(f"Blade Thickness range: {blade_thickness_sorted[0]/self.IN_TO_FT:.4f} to {blade_thickness_sorted[-1]/self.IN_TO_FT:.4f} in")
        print(f"K STAGE SPEED OF SOUND AT INLET  {self.CI:9.2f}")


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
        self.DTDR = np.zeros(shape_mx) # dTheta/dR
        self.DTDZ = np.zeros(shape_mx) # dTheta/dZ
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

    # Removed _spline_deriv_at_points - use spline object directly
    # Removed _spline_interpolate - use spline object directly

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

    def _get_thickness(self, z_coord):
        """
        Interpolates blade thickness at a given Z coordinate using the precomputed spline.
        Replaces _linint for thickness lookup. Handles bounds.

        Args:
            z_coord (float): Z coordinate (in feet) for interpolation.

        Returns:
            float: Interpolated thickness value (in feet), or 0 if out of bounds.
        """
        # Uses the precomputed 1D CubicSpline for thickness
        # Need to handle potential extrapolation if z_coord is outside spline range
        try:
            # Get thickness, fill with 0 if outside bounds
            thickness = self.blade_thickness_spline(z_coord, extrapolate=False)
            # If extrapolate=False returns nan for out of bounds, convert to 0
            return float(np.nan_to_num(thickness))
        except ValueError: # Should not happen with extrapolate=False, but just in case
             return 0.0


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

        if wtfl_calc <= 1e-9: # Avoid division by zero or near-zero
             # If flow is zero, maybe increase velocity slightly? Needs tuning.
             new_wa = wa_current * 1.1 if wa_current > 1e-3 else 1.0 # Start from 1 if current is tiny
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
             print(f"Warning: Large change in WA for orthogonal {i_stream}. "
                   f"Old={wa_current:.4f}, New={new_wa:.4f}, WTFL={wtfl_calc:.4f}")
             # Consider adding logic from original CONTIN for bracketing if needed

        # Original code had error exit (IND=6)
        # if some_divergence_condition:
        #    print(f"ERROR: Convergence failed for streamline {i_stream}, Max WT = {wtfl_calc}")
        #    return None, 6 # Signal failure

        return new_wa, new_ind


    def run_solver(self):
        """Performs the main iterative calculation loop."""

        initial_iter_setting = self.ITER # Store initial iteration count setting
        run_to_convergence = (initial_iter_setting <= 0)

        while self.ITER >= 0: # Loop continues as long as ITER is non-negative
            if run_to_convergence: # Only print headers if running till convergence
                 print(f"\nIteration No. {self.itno:3d}")
                 # Only print detailed headers every few iterations to reduce output
                 if self.itno == 1 or self.itno % 10 == 0:
                     print(f"{' ':6}{'AL':>9}{'RC':>9}{'SM':>9}{'BETA':>9}{'TT':>9}"
                           f"{'SA':>9}{'SB':>9}{'SC':>9}{'SD':>9}")

            self.errori = self.error # Store error from previous iteration
            self.error = 0.0        # Reset max error for current iteration

            # --- Calculate Geometric Parameters and Coefficients (Loops 160-230) ---
            for k in range(self.KMX): # Loop over streamlines
                # Calculate SM (Meridional distance along streamline)
                self.SM[0, k] = 0.0
                for i in range(1, self.MX):
                     ds = math.sqrt((self.Z[i, k] - self.Z[i-1, k])**2 +
                                    (self.R[i, k] - self.R[i-1, k])**2)
                     self.SM[i, k] = self.SM[i-1, k] + ds

                # Fit spline along the streamline (using SM as independent variable)
                # Need Z(SM) and R(SM) splines to get derivatives dZ/dSM, dR/dSM
                if self.SM[-1, k] > 1e-9: # Ensure SM has non-zero length
                    # Check for duplicate SM values which break CubicSpline
                    unique_sm, unique_idx = np.unique(self.SM[:, k], return_index=True)
                    if len(unique_sm) < 2: # Need at least 2 unique points for spline
                         self.AL[:, k].fill(0.0)
                         self.CURV[:, k].fill(0.0)
                    else:
                        spline_z = CubicSpline(unique_sm, self.Z[unique_idx, k], bc_type='natural')
                        spline_r = CubicSpline(unique_sm, self.R[unique_idx, k], bc_type='natural')
                        # Evaluate derivatives at original SM points
                        dZdSM = spline_z.derivative(nu=1)(self.SM[:, k])
                        dRdSM = spline_r.derivative(nu=1)(self.SM[:, k])
                        d2ZdSM2 = spline_z.derivative(nu=2)(self.SM[:, k])
                        d2RdSM2 = spline_r.derivative(nu=2)(self.SM[:, k])

                        # Meridional angle AL = atan(dR/dZ) = atan( (dR/dSM) / (dZ/dSM) )
                        self.AL[:, k] = np.arctan2(dRdSM, dZdSM) # More robust than atan(dR/dZ)

                        # Curvature = (z'r'' - z''r') / (z'^2 + r'^2)^(3/2) where ' = d/dSM
                        denom = (dZdSM**2 + dRdSM**2)**1.5
                        denom[denom < 1e-12] = 1e-12 # Avoid division by zero more strictly
                        self.CURV[:, k] = (dZdSM * d2RdSM2 - d2ZdSM2 * dRdSM) / denom
                else: # Handle zero-length streamline (e.g., first point)
                    self.AL[:, k].fill(0.0) # Or estimate from geometry
                    self.CURV[:, k].fill(0.0) # Curvature is ill-defined

                self.CAL[:, k] = np.cos(self.AL[:, k])
                self.SAL[:, k] = np.sin(self.AL[:, k])

                # Calculate blade angle derivatives (dTheta/dZ and dTheta/dR)
                # DTDZ = dTheta/dZ from the blade angle spline
                self.DTDZ = self.blade_angle_spline.derivative(nu=1)(self.Z[:, k])

                # DTDR = dTheta/dR based on formula involving RB and BETIN (radial effect)
                rinlet = (self.RS[0] + self.RH[0]) / 2.0 # Inlet radius (approx)
                cef_term = math.tan(self.BETIN) / rinlet / (rinlet - self.RB)**2 if abs(rinlet - self.RB) > 1e-9 and rinlet > 0 else 0.0

                for i in range(self.MX):
                    if self.R[i, k] <= self.RB:
                        self.DTDR[i] = 0.0
                    else:
                        self.DTDR[i] = cef_term * (self.R[i, k] - self.RB)**2

                    # Interpolate blade thickness T using the thickness spline
                    t_thick = self._get_thickness(self.Z[i, k]) # Pass Z in feet

                    # Calculate TT (Tangential thickness) and BETA (Relative flow angle)
                    tq = self.R[i, k] * self.DTDR[i] # R * dTheta/dR
                    tp = self.R[i, k] * self.DTDZ[i] # R * dTheta/dZ (Theta from spline)

                    # TT(I,K) = T*SQRT(1.0+TP*TP) - Tangential thickness
                    self.TT[i, k] = t_thick * math.sqrt(1.0 + tp**2)

                    # BETA(I,K) = ATAN(TP*CAL(I,K)+TQ*SAL(I,K)) - Relative flow angle (radians)
                    self.BETA[i, k] = math.atan(tp * self.CAL[i, k] + tq * self.SAL[i, k])
                    self.SBETA[i, k] = math.sin(self.BETA[i, k])
                    self.CBETA[i, k] = math.cos(self.BETA[i, k])

                    # Calculate coefficients SA, SC (depend on geometry and flow angle)
                    self.SA[i, k] = (self.CBETA[i, k]**2 * self.CAL[i, k] * self.CURV[i, k]
                                     - (self.SBETA[i, k]**2 / self.R[i, k] if self.R[i, k] > 1e-9 else 0.0)
                                     + self.SAL[i, k] * self.CBETA[i, k] * self.SBETA[i, k] * self.DTDR[i])

                    self.SC[i, k] = (-self.SAL[i, k] * self.CBETA[i, k]**2 * self.CURV[i, k]
                                     + self.SAL[i, k] * self.CBETA[i, k] * self.SBETA[i, k] * self.DTDZ[i])

                # Calculate velocity derivatives using splines (DWMDM, DWTDM)
                self.AB[:] = self.WA[:, k] * self.CBETA[:, k] # Meridional velocity component?
                self.AC[:] = self.WA[:, k] * self.SBETA[:, k] # Tangential velocity component relative?

                if self.SM[-1, k] > 1e-9: # Check if streamline has length
                    unique_sm, unique_idx = np.unique(self.SM[:, k], return_index=True)
                    if len(unique_sm) < 2:
                         self.DWMDM[:].fill(0.0)
                         self.DWTDM[:].fill(0.0)
                    else:
                        spline_AB = CubicSpline(unique_sm, self.AB[unique_idx], bc_type='natural')
                        spline_AC = CubicSpline(unique_sm, self.AC[unique_idx], bc_type='natural')
                        self.DWMDM[:] = spline_AB.derivative(nu=1)(self.SM[:, k]) # d(W_meridional)/d(SM)
                        self.DWTDM[:] = spline_AC.derivative(nu=1)(self.SM[:, k]) # d(W_tangential_rel)/d(SM)
                else:
                    self.DWMDM[:].fill(0.0)
                    self.DWTDM[:].fill(0.0)

                # Calculate coefficients SB, SD (depend on velocity derivatives)
                for i in range(self.MX):
                    self.SB[i, k] = (self.SAL[i, k] * self.CBETA[i, k] * self.DWMDM[i]
                                     - 2.0 * self.W * self.SBETA[i, k]
                                     + self.DTDR[i] * self.R[i, k] * self.CBETA[i, k] *
                                       (self.DWTDM[i] + 2.0 * self.W * self.SAL[i, k]))

                    self.SD[i, k] = (self.CAL[i, k] * self.CBETA[i, k] * self.DWMDM[i]
                                     + self.DTDZ[i] * self.R[i, k] * self.CBETA[i, k] *
                                       (self.DWTDM[i] + 2.0 * self.W * self.SAL[i, k]))

                    # --- Print intermediate results (if running to convergence) ---
                    if run_to_convergence and k % self.NPRT == 0:
                        # Reduce frequency of printing points along streamline
                        if i == 0: # Print streamline header only once
                             print(f"  Streamline{k+1:3d}")
                        if i % (self.MX // 5 + 1) == 0: # Print roughly 5 points per streamline
                            a_deg = self.AL[i, k] * self.RAD_TO_DEG
                            b_sm_in = self.SM[i, k] / self.IN_TO_FT
                            e_tt_in = self.TT[i, k] / self.IN_TO_FT
                            g_beta_deg = self.BETA[i, k] * self.RAD_TO_DEG
                            print(f"  i={i:3d} "
                                  f"{a_deg:9.3f}{self.CURV[i, k]:9.3f}{b_sm_in:9.3f}"
                                  f"{g_beta_deg:9.3f}{e_tt_in:9.3f}"
                                  f"{self.SA[i, k]:9.3f}{self.SB[i, k]:9.3f}"
                                  f"{self.SC[i, k]:9.3f}{self.SD[i, k]:9.3f}")

            # --- Calculate Blade Surface Velocities (WTR) (Loop 250) ---
            # Calculate every iteration if running to convergence
            if run_to_convergence:
                for k in range(self.KMX):
                    if self.SM[-1, k] > 1e-9: # Check streamline length
                        unique_sm, unique_idx = np.unique(self.SM[:, k], return_index=True)
                        if len(unique_sm) < 2:
                             self.WTR[:, k].fill(0.0)
                             continue

                        # DELBTA = d(TT)/d(SM)
                        spline_TT = CubicSpline(unique_sm, self.TT[unique_idx, k], bc_type='natural')
                        self.DELBTA[:] = spline_TT.derivative(nu=1)(self.SM[:, k])

                        # AB = (R*W + WA*SBETA) * (Circumference / N - TT)
                        # Need DRDM = d(AB)/d(SM)
                        circum_factor = 2.0 * self.PI / self.XN
                        self.AB[:] = ((self.R[:, k] * self.W + self.WA[:, k] * self.SBETA[:, k]) *
                                      (circum_factor * self.R[:, k] - self.TT[:, k]))
                        spline_AB_surf = CubicSpline(unique_sm, self.AB[unique_idx], bc_type='natural')
                        self.DRDM[:] = spline_AB_surf.derivative(nu=1)(self.SM[:, k]) # Store in DRDM

                        # Handle splitter factor SFACT > 1.0
                        ad_temp = np.zeros_like(self.DRDM) # Temporary array like AD in Fortran
                        if self.SFACT > 1.0:
                            circum_factor_s = 2.0 * self.PI / (self.SFACT * self.XN)
                            self.AB[:] = ((self.R[:, k] * self.W + self.WA[:, k] * self.SBETA[:, k]) *
                                          (circum_factor_s * self.R[:, k] - self.TT[:, k]))
                            spline_AB_surf_s = CubicSpline(unique_sm, self.AB[unique_idx], bc_type='natural')
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
                max_contin_attempts = 30 # Increase attempts slightly

                while ind != 0 and contin_attempts < max_contin_attempts:
                    contin_attempts += 1
                    # Update WA along the quasi-orthogonal (k=1 to KMX-1) based on WA[i,0] and coeffs
                    for k in range(1, self.KMX):
                        j = k - 1
                        hr = self.R[i, k] - self.R[i, j]
                        hz = self.Z[i, k] - self.Z[i, j]

                        # Predictor step (WAS)
                        was = self.WA[i, j] * (1.0 + self.SA[i, j] * hr + self.SC[i, j] * hz)

                        # Corrector step uses predicted WAS to estimate influence at k
                        wass = (self.WA[i, j] +
                                was * (self.SA[i, k] * hr + self.SC[i, k] * hz) +
                                self.SB[i, k] * hr + self.SD[i, k] * hz)

                        # Average predictor and corrector
                        self.WA[i, k] = (was + wass) / 2.0

                        # Ensure WA doesn't become negative (physical constraint)
                        if self.WA[i, k] < 0:
                             self.WA[i, k] = 1e-3


                    # Calculate density, pressure, and integrand for weight flow (Loop 340)
                    self.AD_kmx.fill(0.0) # Reset integrand array (size KMX)
                    negative_tip_occured = False
                    for k in range(self.KMX):
                        # Thermodynamics: Calculate Temperature ratio TIP
                        denom_tip = 2.0 * self.CP * self.TEMP
                        if abs(denom_tip) < 1e-9:
                            tip = 1.0
                        else:
                            tip = 1.0 - (self.WA[i, k]**2 + 2.0 * self.W * self.ALM - (self.W * self.R[i, k])**2) / denom_tip

                        # Check for non-physical TIP (negative absolute temperature)
                        if tip < 0.0:
                            # print(f"Warning: Negative TIP encountered at i={i}, k={k}. Halving WA[{i},0].")
                            wa_i0_current *= 0.5
                            self.WA[i, 0] = max(1e-3, wa_i0_current) # Adjust hub velocity, ensure positive
                            ind = 1 # Force re-calculation of WA along orthogonal
                            negative_tip_occured = True
                            break # Exit k loop, restart inner while loop with new WA[i,0]

                        # Calculate TPPIP (related to total conditions?)
                        tppip = 1.0 - (2.0 * self.W * self.ALM - (self.W * self.R[i, k])**2) / denom_tip if abs(denom_tip) > 1e-9 else 1.0

                        # Calculate density with loss term
                        sm_mx_k = self.SM[self.MX-1, k] if self.MX > 0 and self.SM[self.MX-1, k] > 1e-9 else 1.0 # Avoid division by zero
                        loss_term = 0.0
                        if abs(self.AR * tppip * self.TEMP) > 1e-9 and tppip > 1e-9:
                            loss_term = ((tip / tppip)**self.EXPON * self.PLOSS /
                                         (self.AR * tppip * self.TEMP) * 32.17 * # Assuming 32.17 is g
                                         self.SM[i, k] / sm_mx_k)

                        density = tip**self.EXPON * self.RHO - loss_term
                        if density < 0:
                            density = 1e-4

                        # Calculate pressure (PRS)
                        self.PRS[i, k] = density * self.AR * tip * self.TEMP / 32.17 / self.PSI_TO_PSF # Convert psf to psi

                        # Calculate integrand for weight flow (AD(K))
                        # PSI = angle of the quasi-orthogonal (normal to meridional direction)
                        # Use derivatives of hub/shroud curves if available, or approximate
                        # Approximate: Angle of line connecting hub and shroud points
                        dz_qo = self.ZS[i] - self.ZH[i]
                        dr_qo = self.RS[i] - self.RH[i]
                        psi = math.atan2(dr_qo, dz_qo) - self.PI / 2.0 # Angle normal to hub-shroud line

                        # WTHRU = WA(I,K)*CBETA(I,K)*COS(PSI-AL(I,K)) (Velocity normal to orthogonal)
                        wthru = self.WA[i, k] * self.CBETA[i, k] * math.cos(psi - self.AL[i, k])

                        # Effective number of blades (consider splitter)
                        a_blades = self.SFACT * self.XN if self.Z[i, k] < self.ZSPLIT else self.XN

                        # C = passage width = 2*pi*R - N*TT
                        c_width = 2.0 * self.PI * self.R[i, k] - a_blades * self.TT[i, k]
                        if c_width < 0:
                            c_width = 1e-4 # Avoid negative flow

                        # Integrand AD(K) = density * W_normal * Area_element
                        # Area element is C * d(DN)
                        self.AD_kmx[k] = density * wthru * c_width

                    if negative_tip_occured:
                        continue # Restart the inner while loop for this orthogonal 'i'

                    # Integrate AD vs DN to get weight flow
                    self.WTFL = self._cumulative_integral(self.DN[i, :], self.AD_kmx)
                    wtfl_calc = self.WTFL[-1] # Total calculated weight flow for this orthogonal

                    # Check convergence and adjust WA[i, 0] using CONTIN logic
                    if abs(self.WT - wtfl_calc) <= self.WTOLER:
                        ind = 0 # Converged for this orthogonal
                    else:
                        # Call the CONTIN equivalent
                        wa_i0_current, ind = self._contin(wa_i0_current, wtfl_calc, ind, i)
                        if ind == 6: # Check if CONTIN signalled failure
                            print(f"Error: CONTIN failed for orthogonal {i}")
                            self.ITER = -1 # Force exit
                            break # Exit inner while loop
                        self.WA[i, 0] = wa_i0_current # Update hub velocity for next attempt

                    if ind != 0 and contin_attempts >= max_contin_attempts:
                         print(f"Warning: Max CONTIN attempts reached for orthogonal {i}. "
                               f"Calculated WTFL={wtfl_calc:.4f}, Target WT={self.WT:.4f}")
                         ind = 0 # Stop trying for this orthogonal

                if self.ITER < 0: break # Exit outer loop if CONTIN failed

                # --- Update Streamline Positions (DN) ---
                if self.ITER >= 0: # Only update if not exiting due to error
                    # Interpolate DN based on target weight flow distribution BA
                    # Ensure WTFL is monotonically increasing for interpolation
                    if np.all(np.diff(self.WTFL) >= -1e-9): # Allow for small numerical noise
                         # Use linear interpolation as spline might be unstable if WTFL is noisy
                         dn_new = np.interp(self.BA, self.WTFL, self.DN[i, :])
                    else:
                         # Handle non-monotonic WTFL - maybe keep old DN?
                         print(f"Warning: Non-monotonic WTFL encountered for orthogonal {i}. WTFL={self.WTFL}. Skipping DN update.")
                         dn_new = self.DN[i, :].copy() # Keep old values

                    # Calculate change and update DN with correction factor
                    delta = np.abs(dn_new - self.DN[i, :])
                    self.DN[i, :] = (1.0 - self.CORFAC) * self.DN[i, :] + self.CORFAC * dn_new

                    # Update maximum error
                    max_delta_i = np.max(delta) if delta.size > 0 else 0.0
                    if max_delta_i > self.error:
                        self.error = max_delta_i

            # End of loop for quasi-orthogonals (i)
            if self.ITER < 0: break # Exit outer loop if CONTIN failed

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
                print(f"\nConverged after {self.itno} iterations.")
                self.ITER = -1 # Stop iteration
            elif self.error >= self.errori and self.itno > 5: # Check if error increased (allow few initial fluctuations)
                print(f"\nWarning: Error increased or stalled from {self.errori / self.IN_TO_FT:.6f} to {self.error / self.IN_TO_FT:.6f} inches. Stopping.")
                self.ITER = -1 # Stop iteration
            else:
                # Normal iteration decrement (only if initial ITER > 0)
                if not run_to_convergence:
                    self.ITER -= 1

            # Increment iteration number for next loop or final output
            self.itno += 1
            if run_to_convergence:
                 print(f"Iteration {self.itno-1:3d} complete. Max dDN = {self.error / self.IN_TO_FT:10.6f} inches")


            # Safety break if running to convergence and it runs too long
            if run_to_convergence and self.itno > 500: # Arbitrary limit
                 print("Warning: Maximum iterations (500) reached without convergence.")
                 self.ITER = -1 # Force stop


        # --- End of Iteration Loop ---
        print("\n--- Final Results ---")
        self.write_output()

        # Final error print
        print(f"\nFinal Iteration No. {self.itno-1:3d}    Max. Streamline Change = {self.error / self.IN_TO_FT:10.6f} (inches)")


    def write_output(self):
        """Prints the final calculated results."""
        print(f"\n{' ':9}{'DN':>19}{'Z':>19}{'R':>19}"
              f"{'WA':>19}{'PRS':>19}{'WTR':>18}{'RC':>18}")

        # Recalculate final geometry parameters (AL, CURV) if needed for output
        # Optional: Add recalculation here if exact final values are needed

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

    def plot_streamline_parameters(self, streamline_index):
        """
        Plots various calculated parameters along a specified streamline.

        Args:
            streamline_index (int): The index (0-based) of the streamline to plot.
        """
        if not (0 <= streamline_index < self.KMX):
            print(f"Error: Streamline index {streamline_index} out of range (0-{self.KMX-1}).")
            return

        k = streamline_index
        # Ensure SM has valid values before plotting
        if np.any(np.isnan(self.SM[:, k])) or self.SM[-1, k] <= 1e-9:
             print(f"Warning: Insufficient or invalid SM data for streamline {k}. Skipping plot.")
             return

        sm_k = self.SM[:, k] / self.IN_TO_FT # Convert SM to inches for plotting consistency

        params_to_plot = {
            'Meridional Angle (deg)': self.AL[:, k] * self.RAD_TO_DEG,
            'Relative Flow Angle (deg)': self.BETA[:, k] * self.RAD_TO_DEG,
            'Meridional Curvature (1/ft)': self.CURV[:, k],
            'Static Pressure (psi)': self.PRS[:, k],
            'Radius (in)': self.R[:, k] / self.IN_TO_FT,
            'Axial Coordinate Z (in)': self.Z[:, k] / self.IN_TO_FT,
            'Tangential Thickness (in)': self.TT[:, k] / self.IN_TO_FT,
            'Absolute Velocity (ft/s)': self.WA[:, k],
            'Blade Surface Rel. Vel. (ft/s)': self.WTR[:, k]
        }

        num_plots = len(params_to_plot)
        ncols = 3
        nrows = math.ceil(num_plots / ncols)
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5 * nrows), squeeze=False) # Ensure axes is 2D
        axes = axes.flatten() # Flatten to easily iterate

        fig.suptitle(f'Parameters along Streamline {k+1} (Index {k})', fontsize=16)

        plot_count = 0
        for i, (label, data) in enumerate(params_to_plot.items()):
            if np.any(np.isnan(data)):
                 print(f"Warning: NaN values found in '{label}' for streamline {k}. Skipping plot.")
                 continue
            ax = axes[plot_count]
            ax.plot(sm_k, data)
            ax.set_title(label)
            ax.set_xlabel('Meridional Distance SM (in)')
            ax.grid(True)
            plot_count += 1

        # Hide any unused subplots
        for j in range(plot_count, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout
        plt.show()

    def plot_geometry(self):
        """Plots the hub, shroud, and calculated streamline geometries."""
        plt.figure(figsize=(10, 8))
        plt.plot(self.ZH / self.IN_TO_FT, self.RH / self.IN_TO_FT, 'k-', linewidth=2, label='Hub')
        plt.plot(self.ZS / self.IN_TO_FT, self.RS / self.IN_TO_FT, 'k-', linewidth=2, label='Shroud')

        # Plot every NPRT streamline, or fewer if too many
        step = max(1, self.NPRT, self.KMX // 10) # Ensure at least 1, plot NPRT or up to 10 lines
        for k in range(0, self.KMX, step):
             # Check for NaN before plotting
             if not (np.any(np.isnan(self.Z[:, k])) or np.any(np.isnan(self.R[:, k]))):
                 plt.plot(self.Z[:, k] / self.IN_TO_FT, self.R[:, k] / self.IN_TO_FT, '--', label=f'Streamline {k+1}' if k==0 or k==self.KMX-1 else None)

        plt.title('Flowpath Geometry')
        plt.xlabel('Axial Coordinate Z (in)')
        plt.ylabel('Radial Coordinate R (in)')
        plt.legend()
        plt.grid(True)
        plt.axis('equal') # Ensure correct aspect ratio
        plt.show()

    def plot_blade_data(self):
        """Plots the generated blade angle and thickness curves."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Get Z coordinates and corresponding angles/thicknesses from splines
        z_plot = self.blade_angle_spline.x # Use the sorted Z coords stored in spline
        angle_plot_deg = self.blade_angle_spline(z_plot) * self.RAD_TO_DEG
        thickness_plot = self.blade_thickness_spline(z_plot) # Thickness in feet

        # Plot Blade Angle vs Z
        axes[0].plot(z_plot / self.IN_TO_FT, angle_plot_deg)
        axes[0].set_title('Generated Blade Angle')
        axes[0].set_xlabel('Axial Coordinate Z (in)')
        axes[0].set_ylabel('Blade Angle (deg)')
        axes[0].grid(True)

        # Plot Blade Thickness vs Z
        axes[1].plot(z_plot / self.IN_TO_FT, thickness_plot / self.IN_TO_FT)
        axes[1].set_title('Generated Blade Thickness')
        axes[1].set_xlabel('Axial Coordinate Z (in)')
        axes[1].set_ylabel('Thickness (in)')
        axes[1].grid(True)

        plt.tight_layout()
        plt.show()

# --- Main execution ---
if __name__ == "__main__":
    # Create solver instance with desired number of QO lines (MX)
    solver = TurbomachinerySolver(mx_points=30) # Example: Use 30 QO lines

    try:
        solver.setup_from_bezier() # Setup using hardcoded Bezier curves
        solver.run_solver()
        # solver.dump_bcd() # Optional: Dump data if needed

        print("\nRun complete.")

        # Plotting section
        try:
            # Plot geometry and blade data
            solver.plot_geometry()
            solver.plot_blade_data()

            # Ask user for streamline to plot
            while True:
                try:
                    idx_str = input(f"Enter streamline index to plot (0 to {solver.KMX-1}, or 'q' to quit plotting): ")
                    if idx_str.lower() == 'q':
                        break
                    idx = int(idx_str)
                    solver.plot_streamline_parameters(idx)
                except ValueError:
                    print("Invalid input. Please enter an integer index or 'q'.")
                except EOFError:
                     print("\nEOF detected, exiting plotting.")
                     break
                except Exception as plot_err:
                    print(f"Error during plotting: {plot_err}")

        except Exception as e:
            print(f"An error occurred during plotting setup: {e}")

    except KeyboardInterrupt:
        print("\nExecution interrupted by user.")
    except Exception as e:
        print(f"\nAn error occurred during processing: {e}")
        import traceback
        traceback.print_exc()