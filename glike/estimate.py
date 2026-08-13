from glike import *

import cma
import numpy

class Search():
  def __init__(self, x0, bounds = None, precision = 0.05):
    self.names = list(x0.keys())
    self.values = x0.copy()
    if bounds is None:
      bounds = [(0, math.inf) for _ in self.names]
    self.bounds = dict(zip(self.names, bounds))
    self.lrs = {name:0.1 for name in self.names}
    self.precision = precision
  
  def get(self):
    return self.values
  
  def set(self, values):
    self.values = values
  
  def limit(self, name):
    limit = self.bounds[name]
    low = limit[0]; high = limit[1]
    if isinstance(low, str):
      low = eval(low, self.values.copy())
    if isinstance(high, str):
      high = eval(high, self.values.copy())
    return low, high
  
  def up(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = min(high, round(value + step, 5))
    return values
  
  def down(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = round(max(low + 1e-5, value - step), 5)
    return values
  
  def faster(self, name):
    self.lrs[name] = min(0.5, self.lrs[name] * 1.5)
  
  def slower(self, name):
    self.lrs[name] = max(self.precision, self.lrs[name] * 0.5)
  
  def cold(self):
    for name in self.names:
      if self.lrs[name] > self.precision:
        return False
    return True

  def encode_parameters(self, x0):
      """
      Convert user parameters into CMA-ES parameters.
  
      Only transforms ordered times:
          t1 < t2 < t3 < t4
      """
  
      z = dict(x0)
  
      if all(k in x0 for k in ("t1", "t2", "t3", "t4")):
  
          z["dt1"] = np.log(x0["t1"])
          z["dt2"] = np.log(x0["t2"] - x0["t1"])
          z["dt3"] = np.log(x0["t3"] - x0["t2"])
          z["dt4"] = np.log(x0["t4"] - x0["t3"])
  
          del z["t1"]
          del z["t2"]
          del z["t3"]
          del z["t4"]
  
      return z

  def decode_parameters(self, z):
      """
      Convert CMA-ES parameters back into demographic parameters.
      """
  
      x = dict(z)
  
      if all(k in z for k in ("dt1", "dt2", "dt3", "dt4")):
  
          dt1 = np.exp(z["dt1"])
          dt2 = np.exp(z["dt2"])
          dt3 = np.exp(z["dt3"])
          dt4 = np.exp(z["dt4"])
  
          x["t1"] = dt1
          x["t2"] = dt1 + dt2
          x["t3"] = dt1 + dt2 + dt3
          x["t4"] = dt1 + dt2 + dt3 + dt4
  
          del x["dt1"]
          del x["dt2"]
          del x["dt3"]
          del x["dt4"]
  
      return x

def maximize(fun, x0, bounds = None, precision = 0.05, epochs = 20, verbose = False):
  # fun: The objective function to be maximized.
  # x0: the dict of initial parameters, such that the initial output would be fun(**x0)
  # bounds: the list of 2-tuples that defines the boundaries
  # precision: a float that defines the (proportional) step size
  # epochs: an integer that defines the maximum number of epochs
  # verbose: True if intermediate results are printed, False if not
  
  search = Search(x0, bounds = bounds, precision = precision)
  names = list(x0.keys())
  
  y0 = fun(**x0)
  print(str(x0) + " " + str(y0), flush = True)
  
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
        print(" ", flush = True)
        print("x_up: " + str(x_up) + " " + str(y_up), flush = True)
        print("x: " + str(x) + " " + str(y), flush = True)
        print("x_down: " + str(x_down) + " " + str(y_down), flush = True)
        print(" ", flush = True)
      
      if (y_up > max(y_down, y)):
        search.set(x_up)
        search.faster(name)
      elif (y_down > max(y_up, y)):
        search.set(x_down)
        search.faster(name)
      else:
        search.slower(name)
    
    x = search.get()
    y = fun(**x)
    xs.append(x); ys.append(y)
    print(str(x) + " " + str(y), flush = True)
    
    if len(ys) >= 5 and sum(ys[-5:-3]) >= sum(ys[-2:]):
      break
  
  idx = ys.index(max(ys))
  x, y = xs[idx], ys[idx]
  return x, y

def maximize_CMA_ES(fun, x0, bounds = None, precision = 0.05, epochs = 5, verbose = False):
  # fun: Objective function to be maximized
  # x0: dictionary of initial params
  # bounds: Tuple list of paramer bounds
  # precision: Initial search radius
  # epochs: Maximum iterations for CMA-ES
  # verbose: boolean flag to print intermediate progress
  # return best_x: dict of fit params
  # return best_y: log likelihood float

  search = Search(x0, bounds = bounds, precision = precision)
  z0 = search.encode_parameters(x0)

  names = list(z0.keys())

  x_init = np.array([z0[k] for k in names], dtype=float)

  # Check boundary conditions, else make -Inf to Inf
  if bounds is None:
    lower = [-np.inf] * len(names)
    upper = [np.inf] * len(names)
  else:
    # Split boundary pairs
    lower = [b[0] for b in bounds]
    upper = [b[1] for b in bounds]
  # inverse minimzation process of CMA-ES --> maximize likelihood
  def objective(x):
    params = search.decode_parameters(dict(zip(names, x)))

    try:
      y = fun(**params)
      if np.isnan(y):
        return np.inf

      return -y # negative likelihood for minimize
    except Exception as e:
        print("FAILED PARAMETERS:")
        print(params)
        print(e)
        raise
  opts = {
    # "bounds": [lower, upper],
    "maxiter": epochs,
    "verbose": -9
  }

  y0 = fun(**x0)
  print(y0)

  es = cma.CMAEvolutionStrategy(
    x_init,
    precision,
    opts
  )

  generation = 0

  while not es.stop():
    X = es.ask()

    Y = [objective(x) for x in X]

    es.tell(X, Y)

    generation += 1

    if verbose:
      best = dict(zip(names, es.result.xbest))

      print(
        f"Generation {generation:3d}"
        f"  Likelihood = {-es.result.fbest:.6f}"
        
      )

      print(best)
  
  if es.result.xbest is None:
      raise RuntimeError(
          "CMA-ES failed to find a valid solution. "
          "All objective evaluations returned invalid values."
      )

  best_x = search.decode_parameters(
      dict(zip(names, es.result.xbest))
  )

  best_y = -es.result.fbest

  return best_x, best_y
