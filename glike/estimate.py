from glike import *

import cma
import numpy as np
import math


# =============================================================================
# Search
# =============================================================================

class Search():

    def __init__(
        self,
        x0,
        bounds=None,
        precision=0.05,
        parameter_scales=None
    ):
        self.names = list(x0.keys())
        self.values = x0.copy()

        if bounds is None:
            bounds = [(0, math.inf) for _ in self.names]

        self.bounds = dict(zip(self.names, bounds))

        self.lrs = {
            name: 0.1
            for name in self.names
        }

        self.precision = precision

        # Parameter scales used by CMA-ES.
        #
        # These are NOT the same as the bounds. They define the numerical
        # scale of each parameter in CMA-ES space.
        #
        # For example:
        #
        #   N_afr = 14,474
        #   scale = 50,000
        #
        # gives:
        #
        #   normalized N_afr = 14,474 / 50,000 = 0.28948
        #
        # A CMA-ES sigma of 0.05 therefore corresponds to a change of
        # approximately 2,500 individuals.
        self.parameter_scales = (
            parameter_scales
            if parameter_scales is not None
            else {}
        )

    # -------------------------------------------------------------------------
    # Basic Search functionality
    # -------------------------------------------------------------------------

    def get(self):
        return self.values

    def set(self, values):
        self.values = values

    def limit(self, name):
        limit = self.bounds[name]

        low = limit[0]
        high = limit[1]

        if isinstance(low, str):
            low = eval(low, self.values.copy())

        if isinstance(high, str):
            high = eval(high, self.values.copy())

        return low, high

    def up(self, name):

        value = self.values[name]
        lr = self.lrs[name]

        low, high = self.limit(name)

        if value < (low + high) / 2:
            step = (value - low) * lr
        else:
            step = (high - value) * lr

        step = max(step, 1e-5)

        values = self.values.copy()

        values[name] = min(
            high,
            round(value + step, 5)
        )

        return values

    def down(self, name):

        value = self.values[name]
        lr = self.lrs[name]

        low, high = self.limit(name)

        if value < (low + high) / 2:
            step = (value - low) * lr
        else:
            step = (high - value) * lr

        step = max(step, 1e-5)

        values = self.values.copy()

        values[name] = round(
            max(low + 1e-5, value - step),
            5
        )

        return values

    def faster(self, name):
        self.lrs[name] = min(
            0.5,
            self.lrs[name] * 1.5
        )

    def slower(self, name):
        self.lrs[name] = max(
            self.precision,
            self.lrs[name] * 0.5
        )

    def cold(self):

        for name in self.names:

            if self.lrs[name] > self.precision:
                return False

        return True

    # -------------------------------------------------------------------------
    # CMA-ES parameter encoding
    # -------------------------------------------------------------------------

    def encode_parameters(self, x0):
        """
        Convert demographic parameters into the parameterization used
        internally by CMA-ES.

        Ordered times are represented as log-transformed intervals:

            dt1 = log(t1)
            dt2 = log(t2 - t1)
            dt3 = log(t3 - t2)
            dt4 = log(t4 - t3)

        This guarantees that decoded times satisfy:

            t1 < t2 < t3 < t4

        Other parameters are left unchanged at this stage.

        Normalization/scaling is performed separately by
        normalize_parameters().
        """

        z = dict(x0)

        if all(
            k in x0
            for k in ("t1", "t2", "t3", "t4")
        ):

            z["dt1"] = np.log(
                x0["t1"]
            )

            z["dt2"] = np.log(
                x0["t2"] - x0["t1"]
            )

            z["dt3"] = np.log(
                x0["t3"] - x0["t2"]
            )

            z["dt4"] = np.log(
                x0["t4"] - x0["t3"]
            )

            del z["t1"]
            del z["t2"]
            del z["t3"]
            del z["t4"]

        return z

    def decode_parameters(self, z):
        """
        Convert CMA-ES encoded parameters back into demographic
        parameters.

        The inverse of encode_parameters().
        """

        x = dict(z)

        if all(
            k in z
            for k in ("dt1", "dt2", "dt3", "dt4")
        ):

            dt1 = np.exp(
                z["dt1"]
            )

            dt2 = np.exp(
                z["dt2"]
            )

            dt3 = np.exp(
                z["dt3"]
            )

            dt4 = np.exp(
                z["dt4"]
            )

            x["t1"] = dt1

            x["t2"] = (
                dt1
                + dt2
            )

            x["t3"] = (
                dt1
                + dt2
                + dt3
            )

            x["t4"] = (
                dt1
                + dt2
                + dt3
                + dt4
            )

            del x["dt1"]
            del x["dt2"]
            del x["dt3"]
            del x["dt4"]

        return x

    # -------------------------------------------------------------------------
    # Parameter scaling
    # -------------------------------------------------------------------------

    def get_scale(self, name):
        """
        Return the characteristic scale used to normalize a parameter.

        Explicit parameter_scales take precedence.

        If no explicit scale was supplied, finite bounds are used as a
        fallback.

        Finally, a scale of 1.0 is used.

        The latter is appropriate for already-normalized quantities such
        as proportions in [0, 1] and for log-transformed time intervals.
        """

        # -------------------------------------------------------------
        # 1. Explicit user-provided scale
        # -------------------------------------------------------------

        if name in self.parameter_scales:

            scale = self.parameter_scales[name]

            if scale <= 0:
                raise ValueError(
                    f"Scale for {name} must be > 0 "
                    f"(got {scale})"
                )

            return float(scale)

        # -------------------------------------------------------------
        # 2. Finite bounds
        # -------------------------------------------------------------

        if name in self.bounds:

            low, high = self.limit(name)

            if (
                np.isfinite(low)
                and np.isfinite(high)
                and high > low
            ):

                return float(
                    high - low
                )

        # -------------------------------------------------------------
        # 3. Default
        # -------------------------------------------------------------

        return 1.0

    def normalize_parameters(self, z):
        """
        Convert encoded demographic parameters into normalized CMA-ES
        coordinates.

        Each parameter is divided by its characteristic scale.

        Example:

            N_afr = 14,474
            scale = 50,000

            normalized N_afr = 0.28948

        CMA-ES therefore operates on parameters with broadly comparable
        numerical magnitudes.
        """

        normalized = {}

        for name, value in z.items():

            scale = self.get_scale(name)

            normalized[name] = (
                float(value) / scale
            )

        return normalized

    def denormalize_parameters(self, z):
        """
        Convert normalized CMA-ES coordinates back into encoded
        demographic parameters.
        """

        denormalized = {}

        for name, value in z.items():

            scale = self.get_scale(name)

            denormalized[name] = (
                float(value) * scale
            )

        return denormalized


# =============================================================================
# Original coordinate-search optimizer
# =============================================================================

def maximize(
    fun,
    x0,
    bounds=None,
    precision=0.05,
    epochs=20,
    verbose=False
):
    """
    Original coordinate-search optimizer.

    This function is retained unchanged in principle. The normalization
    described above is specifically used by maximize_CMA_ES().
    """

    search = Search(
        x0,
        bounds=bounds,
        precision=precision
    )

    names = list(x0.keys())

    y0 = fun(**x0)

    print(
        str(x0) + " " + str(y0),
        flush=True
    )

    xs = []
    ys = []

    for _ in range(epochs):

        for name in names:

            x = search.get()
            y = fun(**x)

            x_up = search.up(name)
            y_up = fun(**x_up)

            x_down = search.down(name)
            y_down = fun(**x_down)

            if verbose:

                print(
                    " ",
                    flush=True
                )

                print(
                    "x_up: "
                    + str(x_up)
                    + " "
                    + str(y_up),
                    flush=True
                )

                print(
                    "x: "
                    + str(x)
                    + " "
                    + str(y),
                    flush=True
                )

                print(
                    "x_down: "
                    + str(x_down)
                    + " "
                    + str(y_down),
                    flush=True
                )

                print(
                    " ",
                    flush=True
                )

            if (
                y_up
                > max(
                    y_down,
                    y
                )
            ):

                search.set(x_up)
                search.faster(name)

            elif (
                y_down
                > max(
                    y_up,
                    y
                )
            ):

                search.set(x_down)
                search.faster(name)

            else:

                search.slower(name)

        x = search.get()
        y = fun(**x)

        xs.append(x)
        ys.append(y)

        print(
            str(x) + " " + str(y),
            flush=True
        )

        if (
            len(ys) >= 5
            and sum(ys[-5:-3])
            >= sum(ys[-2:])
        ):

            break

    idx = ys.index(
        max(ys)
    )

    x = xs[idx]
    y = ys[idx]

    return x, y


# =============================================================================
# Invalid-parameter handler
# =============================================================================

def invalid(reason, params):

    print(
        "INVALID PARAMETERS:"
    )

    print(
        reason
    )

    print(
        params
    )

    return np.inf


# =============================================================================
# CMA-ES log-space parameter transform
# =============================================================================

class CMAParameterTransform:
    """
    Map demographic parameters to a common CMA-ES coordinate system [0, 1].

    The transformation is based on one quantity for every parameter:

        d = x - lower_bound

    and optimizes log(d).

    For finite bounds L < x < U:

        u = [log(x - L) - log(U - L)] / [log(U - L) - log(U - L)]

    More explicitly, since the lower endpoint itself is singular under the
    logarithm, the implementation uses a small positive distance from L as
    the lower log-space endpoint. This gives a numerically stable mapping
    whose interior is strictly inside (0, 1).

    A cleaner equivalent representation used below is to define a positive
    log-space interval [log(delta_min), log(delta_max)] where delta_min is
    chosen from the supplied starting point and numerical precision. CMA-ES
    then operates on the normalized position within that interval.

    For an unbounded upper limit, an explicit parameter scale defines the
    reference distance and the same log-distance coordinate is mapped to
    (0, 1) with a logistic transform:

        u = sigmoid(log((x-L)/scale)).

    Thus all CMA-ES coordinates have the same interpretation regardless of
    the original units of the demographic parameter.
    """

    EPS = 1e-12

    def __init__(self, x0, bounds=None, parameter_scales=None):
        self.x0 = {k: float(v) for k, v in x0.items()}
        self.names = list(x0.keys())

        if bounds is None:
            bounds = [(0, math.inf) for _ in self.names]

        self.bounds = dict(zip(self.names, bounds))
        self.parameter_scales = parameter_scales or {}

        self._validate_bounds()

    def _eval_bound(self, value, params):
        if isinstance(value, str):
            return float(eval(value, {"__builtins__": {}}, dict(params)))
        return float(value)

    def _raw_bounds(self, name, params):
        low, high = self.bounds[name]
        return self._eval_bound(low, params), self._eval_bound(high, params)

    def _ordered_lower_bound(self, name, params, low):
        """Enforce strict ordering for consecutively numbered time parameters."""
        if name.startswith("t") and name[1:].isdigit():
            idx = int(name[1:])
            if idx > 1:
                previous = f"t{idx - 1}"
                if previous in params:
                    low = max(low, float(params[previous]))
        return low

    def _validate_bounds(self):
        # Check static bounds and the supplied starting point. Dynamic bounds
        # are checked after earlier parameters have been decoded.
        for name in self.names:
            low, high = self.bounds[name]
            if not isinstance(low, str) and not np.isfinite(float(low)):
                raise ValueError(f"Lower bound for {name} must be finite.")
            if not isinstance(high, str) and float(high) <= float(low):
                raise ValueError(f"Upper bound for {name} must be greater than lower bound.")

    def _scale_for_unbounded(self, name, lower):
        if name in self.parameter_scales:
            scale = float(self.parameter_scales[name])
        else:
            # A scale is only needed when the upper bound is infinite.
            # One unit of distance from the lower bound is the neutral default.
            scale = 1.0

        if scale <= 0 or not np.isfinite(scale):
            raise ValueError(
                f"Scale for {name} must be finite and > 0 (got {scale})."
            )
        return scale

    def _bounds_for_name(self, name, decoded):
        low, high = self._raw_bounds(name, decoded)
        low = self._ordered_lower_bound(name, decoded, low)
        if np.isfinite(high) and high <= low:
            raise ValueError(
                f"Invalid bounds for {name}: lower={low}, upper={high}."
            )
        return low, high

    @staticmethod
    def _logit(u):
        u = min(max(float(u), CMAParameterTransform.EPS), 1.0 - CMAParameterTransform.EPS)
        return math.log(u / (1.0 - u))

    @staticmethod
    def _sigmoid(z):
        z = float(z)
        if z >= 0:
            ez = math.exp(-z)
            return 1.0 / (1.0 + ez)
        ez = math.exp(z)
        return ez / (1.0 + ez)

    def encode(self, x):
        """
        Convert original parameters to CMA-ES coordinates in (0, 1).

        Parameters are processed in their supplied order. This is important
        for ordered times because t2, t3, ... use the previously decoded time
        as their effective lower bound.
        """
        decoded = {}
        u = {}

        for name in self.names:
            value = float(x[name])
            low, high = self._bounds_for_name(name, decoded)
            distance = value - low

            if not distance > 0:
                raise ValueError(
                    f"Initial value for {name} must be strictly above its "
                    f"lower bound: value={value}, lower={low}."
                )

            if np.isfinite(high):
                upper_distance = high - low
                log_d = math.log(distance)
                log_upper = math.log(upper_distance)

                # The lower endpoint of the finite log interval is chosen
                # from machine-scale precision relative to the interval.
                log_lower = math.log(max(upper_distance * self.EPS, np.finfo(float).tiny))

                ui = (log_d - log_lower) / (log_upper - log_lower)
                if not (0.0 < ui < 1.0):
                    raise ValueError(
                        f"Initial value for {name} lies outside its log-scaled "
                        f"finite bounds: value={value}, bounds=({low}, {high})."
                    )
            else:
                scale = self._scale_for_unbounded(name, low)
                ui = self._sigmoid(math.log(distance / scale))

            u[name] = ui
            decoded[name] = value

        return u

    def decode(self, u):
        """
        Convert CMA-ES coordinates in (0, 1) back to original parameters.

        No clipping or model-specific repair is performed. Values outside
        (0, 1) are rejected, making bounds hard constraints.
        """
        params = {}

        for name in self.names:
            ui = float(u[name])
            if not (0.0 < ui < 1.0):
                raise ValueError(
                    f"CMA-ES coordinate for {name} must lie strictly within "
                    f"(0, 1); got {ui}."
                )

            low, high = self._bounds_for_name(name, params)

            if np.isfinite(high):
                upper_distance = high - low
                log_lower = math.log(max(upper_distance * self.EPS, np.finfo(float).tiny))
                log_upper = math.log(upper_distance)
                log_distance = log_lower + ui * (log_upper - log_lower)
                value = low + math.exp(log_distance)
            else:
                scale = self._scale_for_unbounded(name, low)
                value = low + scale * math.exp(self._logit(ui))

            params[name] = float(value)

        return params

    def validate(self, params):
        """Return True only when params satisfy the original hard bounds."""
        checked = {}
        for name in self.names:
            low, high = self._bounds_for_name(name, checked)
            value = float(params[name])
            if value <= low:
                return False
            if np.isfinite(high) and value >= high:
                return False
            checked[name] = value
        return True


# =============================================================================
# CMA-ES optimizer
# =============================================================================

def maximize_CMA_ES(
    fun,
    x0,
    bounds=None,
    precision=0.05,
    epochs=20,
    verbose=False,
    model=None,
    parameter_scales=None
):
    """
    Maximize an objective function using CMA-ES in a common log-scaled space.

    CMA-ES operates exclusively on normalized coordinates in (0, 1). Every
    parameter uses the same mathematical transformation: log of its distance
    above its lower bound. Finite upper bounds define the upper end of the
    log-space interval; infinite upper bounds use parameter_scales as the
    reference distance in the same log transform.

    No model-specific clipping or parameter repair is performed. A candidate
    that violates a hard constraint is assigned an infinite minimization
    objective and therefore cannot become the best solution.
    """

    transform = CMAParameterTransform(
        x0=x0,
        bounds=bounds,
        parameter_scales=parameter_scales
    )

    names = list(x0.keys())
    u0_dict = transform.encode(x0)
    u0 = np.array([u0_dict[name] for name in names], dtype=float)

    if verbose:
        print("\nInitial parameters:")
        print(x0)
        print("\nInitial log-scaled CMA-ES coordinates:")
        print(u0_dict)

    def objective(u):
        try:
            u_dict = dict(zip(names, u))

            # Hard normalized-domain constraint. No clipping.
            if any(not (0.0 < float(v) < 1.0) for v in u):
                return np.inf

            params = transform.decode(u_dict)

            if not transform.validate(params):
                return np.inf

            y = float(fun(**params))

            if not np.isfinite(y):
                return np.inf

            return -y

        except Exception as e:
            if verbose:
                print("FAILED PARAMETERS:")
                try:
                    print(transform.decode(dict(zip(names, u))))
                except Exception:
                    print(dict(zip(names, u)))
                print(f"ERROR: {repr(e)}")
            return np.inf

    y0 = float(fun(**x0))
    print(y0)

    # CMA-ES is constrained to the common normalized coordinate system.
    # The explicit objective-domain check above remains the final hard
    # constraint; pycma's bounds option prevents excursions beyond [0, 1].
    opts = {
        "maxiter": epochs,
        "bounds": [0.0, 1.0],
        "verbose": -9
    }

    es = cma.CMAEvolutionStrategy(
        u0,
        precision,
        opts
    )

    generation = 0

    while not es.stop():
        X = es.ask()
        Y = [objective(u) for u in X]
        es.tell(X, Y)
        generation += 1

        if verbose:
            if es.result.xbest is not None and np.isfinite(es.result.fbest):
                best_u = dict(zip(names, es.result.xbest))
                best_x = transform.decode(best_u)
                print(
                    f"Generation {generation:3d}"
                    f"  Likelihood = {-es.result.fbest:.6f}"
                )
                print(best_x)
            else:
                print(
                    f"Generation {generation:3d}"
                    "  No valid CMA-ES solution yet."
                )

    if es.result.xbest is None or not np.isfinite(es.result.fbest):
        raise RuntimeError("CMA-ES failed to find a valid solution.")

    best_u = dict(zip(names, es.result.xbest))
    best_x = transform.decode(best_u)
    best_y = -float(es.result.fbest)

    if verbose:
        print("\nFinal log-scaled CMA-ES coordinates:")
        print(best_u)
        print("\nFinal demographic parameters:")
        print(best_x)
        print(f"\nFinal likelihood = {best_y:.6f}")

    return best_x, best_y
