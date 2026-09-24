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
# CMA-ES optimizer
# =============================================================================

def maximize_CMA_ES(
    fun,
    x0,
    bounds=None,
    precision=0.05,
    epochs=5,
    verbose=False,
    model=None,
    parameter_scales=None
):
    """
    Maximize an objective function using CMA-ES.

    Parameters
    ----------
    fun : callable
        Objective function to maximize.

    x0 : dict
        Initial demographic parameters.

    bounds : list of tuples, optional
        Parameter bounds in the ORIGINAL parameterization.

    precision : float
        Initial CMA-ES standard deviation in NORMALIZED parameter space.

        For example:

            precision = 0.05

        means that the initial CMA-ES search distribution has a standard
        deviation of 0.05 normalized units for every parameter.

    epochs : int
        Maximum number of CMA-ES generations.

    verbose : bool
        Print intermediate progress.

    model : str, optional
        Demographic model name.

    parameter_scales : dict, optional
        Characteristic scale for each parameter.

        Example:

            {
                "N_afr": 50000,
                "N_eur": 50000,
                "N_asia": 50000,
                "r1": 1.0,
                "r2": 1.0,
                "r3": 1.0,
                "gr": 1.0
            }

        Parameters not included in this dictionary fall back to their
        finite bound range or a scale of 1.0.
    """

    # -------------------------------------------------------------------------
    # Construct Search object
    # -------------------------------------------------------------------------

    if parameter_scales is None:
        parameter_scales = {}

        for name in x0:
            if name.startswith("N_"):
                parameter_scales[name] = 50000
            elif name.startswith("r"):
                parameter_scales[name] = 1.0
            elif name.startswith("m"):
                parameter_scales[name] = 1.0
            elif name.startswith("gr"):
                parameter_scales[name] = 1.0

    search = Search(
        x0,
        bounds=bounds,
        precision=precision,
        parameter_scales=parameter_scales
    )

    # -------------------------------------------------------------------------
    # Encode original parameters
    # -------------------------------------------------------------------------

    z0 = search.encode_parameters(
        x0
    )

    # -------------------------------------------------------------------------
    # Normalize parameters for CMA-ES
    # -------------------------------------------------------------------------

    q0 = search.normalize_parameters(
        z0
    )

    names = list(
        q0.keys()
    )

    # CMA-ES initial vector
    x_init = np.array(
        [
            q0[k]
            for k in names
        ],
        dtype=float
    )

    # -------------------------------------------------------------------------
    # Print initial parameterization
    # -------------------------------------------------------------------------

    if verbose:

        print(
            "\nInitial parameters:"
        )

        print(
            x0
        )

        print(
            "\nEncoded parameters:"
        )

        print(
            z0
        )

        print(
            "\nNormalized CMA-ES parameters:"
        )

        print(
            q0
        )

        print(
            "\nCMA-ES parameter scales:"
        )

        for name in names:

            print(
                f"  {name}: "
                f"{search.get_scale(name)}"
            )

        print(
            ""
        )

    # -------------------------------------------------------------------------
    # Objective function
    # -------------------------------------------------------------------------

    def objective(q):
        """
        Convert normalized CMA-ES coordinates back into the original
        demographic parameters and evaluate the likelihood.
        """

        # -------------------------------------------------------------
        # Normalized CMA-ES coordinates
        # -------------------------------------------------------------

        q_dict = dict(
            zip(
                names,
                q
            )
        )

        # -------------------------------------------------------------
        # Undo parameter scaling
        # -------------------------------------------------------------

        z = search.denormalize_parameters(
            q_dict
        )

        # -------------------------------------------------------------
        # Decode ordered times
        # -------------------------------------------------------------

        params = search.decode_parameters(
            z
        )

        # -------------------------------------------------------------
        # Convert numpy scalars to Python floats
        # -------------------------------------------------------------

        params = {
            k: float(v)
            for k, v in params.items()
        }

        # -------------------------------------------------------------
        # Parameter validity
        # -------------------------------------------------------------

        if (
            model == "3G09"
            or model == "3G09_no_m"
        ):

            if params["N_anc"] <= 0:
                params["N_anc"] = 1e-5

            if params["N_yri"] <= 0:
                params["N_yri"] = 1e-5

            if params["N_ooa"] <= 0:
                params["N_ooa"] = 1e-5

            if params["N_ceu"] <= 0:
                params["N_ceu"] = 1e-5

            if params["N_chb"] <= 0:
                params["N_chb"] = 1e-5

            if params["gr_ceu"] < 0:
                params["gr_ceu"] = 1e-5

            if params["gr_ceu"] > 1:
                params["gr_ceu"] = 1 - 1e-5

            if params["gr_chb"] < 0:
                params["gr_chb"] = 1e-5

            if params["gr_chb"] > 1:
                params["gr_chb"] = 1 - 1e-5

            if params["t1"] <= 0:
                params["t1"] = 1e-5

            if params["t2"] <= params["t1"]:
                params["t2"] = params["t1"] + 1e-5

            if params["t3"] <= params["t2"]:
                params["t3"] = params["t2"] + 1e-5

        elif model == "3I21":

            if params["N_yri"] <= 0:
                params["N_yri"] = 1e-5

            if params["N_ceu"] <= 0:
                params["N_ceu"] = 1e-5

            if params["N_nea"] <= 0:
                params["N_nea"] = 1e-5

            if params["m1"] < 0:
                params["m1"] = 1e-5

            if params["m1"] > 1:
                params["m1"] = 1 - 1e-5

            if params["t1"] <= 0:
                params["t1"] = 1e-5

            if params["t2"] <= params["t1"]:
                params["t2"] = params["t1"] + 1e-5

            if params["t3"] <= params["t2"]:
                params["t3"] = params["t2"] + 1e-5

            if params["t4"] <= params["t3"]:
                params["t4"] = params["t3"] + 1e-5

        elif model == "4A21":

            if params["t1"] <= 0:
                params["t1"] = 1e-5

            if params["t2"] <= params["t1"]:
                params["t2"] = params["t1"] + 1e-5

            if params["t3"] <= params["t2"]:
                params["t3"] = params["t2"] + 1e-5

            if params["t4"] <= params["t3"]:
                params["t4"] = params["t3"] + 1e-5

            if params["t5"] <= params["t4"]:
                params["t5"] = params["t4"] + 1e-5

            if params["t6"] <= params["t5"]:
                params["t6"] = params["t5"] + 1e-5

            if params["r1"] <= 0:
                params["r1"] = 1e-5

            if params["r2"] <= 0:
                params["r2"] = 1e-5

            if params["r3"] <= 0:
                params["r3"] = 1e-5

            if params["N_ana"] <= 0:
                params["N_ana"] = 1e-5

            if params["N_neo"] <= 0:
                params["N_neo"] = 1e-5

            if params["N_whg"] <= 0:
                params["N_whg"] = 1e-5

            if params["N_bronze"] <= 0:
                params["N_bronze"] = 1e-5

            if params["N_yam"] <= 0:
                params["N_yam"] = 1e-5

            if params["N_ehg"] <= 0:
                params["N_ehg"] = 1e-5

            if params["N_chg"] <= 0:
                params["N_chg"] = 1e-5

            if params["N_ne"] <= 0:
                params["N_ne"] = 1e-5

            if params["N_wa"] <= 0:
                params["N_wa"] = 1e-5

            if params["N_ooa"] <= 0:
                params["N_ooa"] = 1e-5

            if params["gr"] < 0:
                params["gr"] = 1e-5

            if params["gr"] > 1:
                params["gr"] = 1 - 1e-5

        else:

            if params["t1"] <= 0:
                params["t1"] = 1e-5

            if params["t2"] <= params["t1"]:
                params["t2"] = params["t1"] + 1e-5

            if params["t3"] <= params["t2"]:
                params["t3"] = params["t2"] + 1e-5

            if params["t4"] <= params["t3"]:
                params["t4"] = params["t3"] + 1e-5

            if params["r1"] <= 0:
                params["r1"] = 1e-5

            if params["r2"] <= 0:
                params["r2"] = 1e-5

            if params["r3"] <= 0:
                params["r3"] = 1e-5

            if params["N_admix"] <= 0:
                params["N_admix"] = 1e-5

            if params["N_afr"] <= 0:
                params["N_afr"] = 1e-5

            if params["N_eur"] <= 0:
                params["N_eur"] = 1e-5

            if params["N_asia"] <= 0:
                params["N_asia"] = 1e-5

            if params["N_pol"] <= 0:
                params["N_pol"] = 1e-5

            if params["N_aa"] <= 0:
                params["N_aa"] = 1e-5

            if params["N_ooa"] <= 0:
                params["N_ooa"] = 1e-5

            if params["N_anc"] <= 0:
                params["N_anc"] = 1e-5

            if params["gr"] < 0:
                params["gr"] = 1e-5

            if params["gr"] >= 1:
                params["gr"] = 1 - 1e-5

        # -------------------------------------------------------------
        # Evaluate likelihood
        # -------------------------------------------------------------

        try:

            y = fun(
                **params
            )

            if not np.isfinite(y):
                return np.inf

            # CMA-ES minimizes, while we want to maximize.
            return -y

        except Exception as e:

            if verbose:

                print(
                    "FAILED PARAMETERS:"
                )

                print(
                    params
                )

                print(
                    f"ERROR: {repr(e)}"
                )

            return np.inf

    # -------------------------------------------------------------------------
    # CMA-ES options
    # -------------------------------------------------------------------------

    opts = {
        "maxiter": epochs,
        "verbose": -9
    }

    # -------------------------------------------------------------------------
    # Evaluate initial likelihood
    # -------------------------------------------------------------------------

    y0 = fun(
        **x0
    )

    print(
        y0
    )

    # -------------------------------------------------------------------------
    # Run CMA-ES
    #
    # IMPORTANT:
    #
    # x_init is normalized.
    #
    # precision is therefore also expressed in normalized units.
    # -------------------------------------------------------------------------

    es = cma.CMAEvolutionStrategy(
        x_init,
        precision,
        opts
    )

    generation = 0

    while not es.stop():

        # -------------------------------------------------------------
        # Generate candidate solutions in normalized space
        # -------------------------------------------------------------

        X = es.ask()

        # -------------------------------------------------------------
        # Evaluate candidates
        # -------------------------------------------------------------

        Y = [
            objective(x)
            for x in X
        ]

        # -------------------------------------------------------------
        # Update CMA-ES
        # -------------------------------------------------------------

        es.tell(
            X,
            Y
        )

        generation += 1

        # -------------------------------------------------------------
        # Optional progress reporting
        # -------------------------------------------------------------

        if verbose:

            if es.result.xbest is not None:

                best_normalized = dict(
                    zip(
                        names,
                        es.result.xbest
                    )
                )

                best_encoded = (
                    search.denormalize_parameters(
                        best_normalized
                    )
                )

                best = (
                    search.decode_parameters(
                        best_encoded
                    )
                )

                print(
                    f"Generation {generation:3d}"
                    f"  Likelihood = "
                    f"{-es.result.fbest:.6f}"
                )

                print(
                    best
                )

            else:

                print(
                    f"Generation {generation:3d}"
                    "  No valid CMA-ES solution yet."
                )

    # -------------------------------------------------------------------------
    # Check for valid result
    # -------------------------------------------------------------------------

    if es.result.xbest is None:

        raise RuntimeError(
            "CMA-ES failed to find a valid solution."
        )

    # -------------------------------------------------------------------------
    # Convert best CMA-ES solution back to original parameters
    # -------------------------------------------------------------------------

    best_normalized = dict(
        zip(
            names,
            es.result.xbest
        )
    )

    best_encoded = (
        search.denormalize_parameters(
            best_normalized
        )
    )

    best_x = (
        search.decode_parameters(
            best_encoded
        )
    )

    # -------------------------------------------------------------------------
    # Best likelihood
    # -------------------------------------------------------------------------

    best_y = (
        -es.result.fbest
    )

    # -------------------------------------------------------------------------
    # Final verbose output
    # -------------------------------------------------------------------------

    if verbose:

        print(
            "\nFinal normalized parameters:"
        )

        print(
            best_normalized
        )

        print(
            "\nFinal demographic parameters:"
        )

        print(
            best_x
        )

        print(
            f"\nFinal likelihood = "
            f"{best_y:.6f}"
        )

    return best_x, best_y