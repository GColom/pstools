#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

template<typename T> inline

int signum(T x) {return (x > T(0)) - (x <  T(0));}

template<typename T>

auto midpoint(std::vector<T> time_series)
{
  auto mid = time_series.size()/2;
  return time_series.begin() + mid;
}

static int timeprecision =
    8;  // 8 is the RADAR5 output precision for the time variable.
static int PSprecision =
    8;  // 8 is the RADAR5 output precision for phase space variables.

// Class implementing a single point in a dynamical system's Phase Space.
class PSpoint {
 public:
  double t, x, y;

  // Constructors
  constexpr PSpoint(double mt, double mx, double my) : t(mt), x{mx}, y{my} {}
  PSpoint() {}

  // Friend operators declarations
  friend std::ostream& operator<<(std::ostream& o, PSpoint p);
  friend std::istream& operator>>(std::ostream& o, PSpoint p);
  friend bool operator==(const PSpoint&, const PSpoint&);

  PSpoint& operator+=(const PSpoint& q) {
    this->t += q.t;
    this->x += q.x;
    this->y += q.y;
    return *this;
  }

  PSpoint& operator-=(const PSpoint& q) {
    this->t -= q.t;
    this->x -= q.x;
    this->y -= q.y;
    return *this;
  }

  PSpoint& operator/(const double d) {
    if (d == 0) throw(std::runtime_error("Division by zero error."));

    this->t /= d;
    this->x /= d;
    this->y /= d;

    return *this;
  }
};

// Insertion/Extraction operators
std::ostream& operator<<(std::ostream& o, PSpoint p) {
  o << std::setprecision(timeprecision) << p.t << '\t';
  o << std::setprecision(PSprecision) << p.x << '\t' << p.y << std::endl;
  return o;
}

std::istream& operator>>(std::istream& i, PSpoint& p) {
  double _t, _x, _y;
  i >> _t;
  i >> _x;
  i >> _y;
  p = PSpoint{_t, _x, _y};
  return i;
}

// Lesser operator, for sorting
bool operator<(const PSpoint& p, const PSpoint& q) { return p.t < q.t; }
// Lesser operator with const double, for usage in std::lower_bound
bool operator<(const PSpoint& p, const double t) { return p.t < t; }
//== operator, for equality check
bool operator==(const PSpoint& p, const PSpoint& q) {
  return (p.t == q.t) and (p.x == q.x) and (p.y == q.y);
}
//+ operator, for summing points
PSpoint operator+(const PSpoint& p, const PSpoint& q) {
  return PSpoint(p.t + q.t, p.x + q.x, p.y + q.y);
}
//- operator, for point subtraction
PSpoint operator-(const PSpoint& p, const PSpoint& q) {
  return PSpoint(p.t - q.t, p.x - q.x, p.y - q.y);
}

// Class implementing the full trajectory of a 2-dimensional dynamical system
class trajectory {
  std::vector<PSpoint> points;

  // Deduce correct type for vector of iterators marking the limits of
  // periods in the trajectory.
  typedef decltype(points.begin()) traj_iterator;

 public:
  trajectory(char* path) {
    std::ifstream filein(path);
    if (filein.good()) {
      filein >> *this;
    } else {
      throw std::runtime_error(
          "Bad ifstream! You possibly tried to open a non existent record.");
    }
  }
  trajectory() {}

  trajectory(std::vector<PSpoint> pts) { points = pts; }

  auto begin() { return points.begin(); }

  auto end() { return points.end(); }

  void push_back(double t, double x, double y) {
    PSpoint p{t, x, y};
    points.push_back(p);
  }

  void add_point(double t, double x, double y) {
    PSpoint p{t, x, y};
    this->push_back(p);
    std::sort(points.begin(), points.end());
  }

  void push_back(PSpoint p) { points.push_back(p); }

  void add_point(PSpoint p) {
    push_back(p);
    std::sort(points.begin(), points.end());
  }

  void sort() { std::sort(points.begin(), points.end()); }

  void pop_back() { points.pop_back(); }

  friend std::ostream& operator<<(std::ostream& o, trajectory t);
  friend std::istream& operator>>(std::istream& i, trajectory& t);
  friend bool operator==(const trajectory&, const trajectory&);

  const PSpoint& operator[](const int& i) const { return points[i]; }

  PSpoint& operator[](const int& i) { return points[i]; }

  void clear() { this->points.clear(); }

  int size() const { return points.size(); }

  double t_end() const {
    auto maxiter = std::max_element(points.begin(), points.end());
    return maxiter->t;
  }

  std::vector<traj_iterator> period_partition(double tau)
  // Return a vector of iterators, so that the part of trajectory
  // between any two consecutive ones corresponds to a time tau.
  {
    std::vector<traj_iterator> res;
    // this->begin() always marks the beginning of the first period.
    res.push_back(this->begin());

    int n_periods = this->t_end() / tau;

    std::vector<double> periods;

    for (int i = 1; i <= n_periods; ++i) periods.push_back(i * tau);

    for (auto p : periods)
      res.push_back(std::lower_bound(this->begin(), this->end(), p));

    return res;
  }

  std::vector<traj_iterator> loop_partition(double x)
  // Return a vector of iterators, so that the part of trajectory
  // between any two consecutive ones corresponds to a full revolution.
  {
    std::vector<traj_iterator> res;
    // this->begin() always marks the first record after the first crossing.
    res.push_back(this->begin());
    traj_iterator it = this->begin();

    while (it != this->end()) {
      it = std::adjacent_find(it, this->end(), [&x](PSpoint cur, PSpoint next) {
        return (cur.x < x) && (next.x >= x);
      });
      if (it != this->end()) {
        // At this point the it marks the first of the two PSpoints where
        // crossing occurs. We must increment it so that it points the first
        // element after the crossing. This way any ranging we might do between
        // any two elements of res will yield the portion of orbit strictly
        // after the first and before the last iterator provided.
        ++it;
        res.push_back(it);
      } else
        break;
    }

    return res;
  }

  double shoelace_area(traj_iterator beg, traj_iterator end) {
    std::vector<double> summands(end - beg - 1, 0.0);
    std::transform(beg, end - 1, beg + 1, summands.begin(),
                   [](PSpoint cur, PSpoint next) {
                     return (cur.y + next.y) * (cur.x - next.x);
                   });
    return std::abs(0.5 *
                    std::accumulate(summands.begin(), summands.end(), 0.0));
  }

  std::vector<double> period_areas(double tau)
  // Return the area enclosed within the trajectory described for each period.
  {
    auto period_markers = this->period_partition(tau);

    int n_markers = period_markers.size();

    std::vector<double> areas(n_markers - 1, 0.0);

    for (int i = 0; i < n_markers - 1; ++i)
      areas[i] = shoelace_area(period_markers[i], period_markers[i + 1]);
    return areas;
  }

  std::vector<double> loop_areas(double x)
  // Return the area enclosed within each revolution.
  // A single -1 in the vector indicates that the system has never crossed the
  // section This can be most likely due to some kind of confinement or
  // overdamping.
  {
    auto loop_markers = this->loop_partition(x);

    int n_markers = loop_markers.size();

    std::vector<double> areas(n_markers - 1, 0.0);

    for (int i = 0; i < n_markers - 1; ++i)
      areas[i] = shoelace_area(loop_markers[i], loop_markers[i + 1]);

    if (areas.size() == 0) areas.push_back(-1.);

    return areas;
  }

  double loop_area_average(double x, double discard_percentage)
  // Return the average of revolution areas, calculated over the last
  // T.size()*(1-discard_percentage revolutions)
  //
  // A single -1 in the vector indicates that the system has never crossed the
  // section This can be most likely due to some kind of confinement or
  // overdamping.
  {
    std::vector<double> areas = this->loop_areas(x);

    int n_areas = areas.size();

    if ((n_areas == 1) && (areas[0] == -1.)) return -1.;

    int head = int(std::round(discard_percentage * n_areas));

    auto beg = areas.begin() + head;

    double n_avg = std::distance(beg, areas.end());

    return std::accumulate(beg, areas.end(), 0.0) / n_avg;
  }


  double loop_area_asymptote(double x, double discard_percentage)
  // Return the asymptotic estimate of the area per loop, discarding a
  // percentage of the initial areas. The estimate is carried out via two
  // successive analytical least squares fits.
  // 
  // The fit is performed onto a function f(x) = A - B exp(-Cx), the return
  // of the function is a vector containing A, B and C in this order.
  //
  // A single -1 in the vector indicates that the system has never crossed the
  // section This can be most likely due to some kind of confinement or
  // overdamping.
  {
    std::vector<double> y = this->loop_areas(x);

    int n_areas = y.size();

    if ((n_areas == 1) && (y[0] == -1.)) return -1.;

    int offset = int(std::round(discard_percentage * n_areas));

    auto beg = y.begin() + offset;

    double n_records = std::distance(beg, y.end());

    std::vector<double> log_deltas(n_records-1, 0.);

    int trend = signum(*y.end()-*midpoint(y));

    std::cout<<trend<<std::endl;

    std::adjacent_difference(beg, y.end(), log_deltas.begin(), [](double a, double b){return std::log(std::abs(b-a));});

    double x_square_avg   = (n_areas*(n_areas+1)*(2*n_areas+1) - offset*(offset+1)*(2*offset+1))/(6*n_records);
    double x_avg          = (n_areas*(n_areas+1) - offset*(offset+1))/(2*n_records);
    double log_deltas_avg = std::accumulate(log_deltas.begin(), log_deltas.end(), 0.0)/(n_records-1);
    double x_log_delta_cov     = 0.0;
    
    int x_it = offset;
    for(auto it_log:log_deltas)
    {
      x_log_delta_cov += it_log * x_it / n_records;
      x_it++;
    }

    double C = - (x_log_delta_cov - x_avg*log_deltas_avg)/(x_square_avg - x_avg*x_avg);
    double log_D = (x_log_delta_cov/x_avg) + (x_square_avg * C)/x_avg;

    double D = std::exp(log_D);
    double B = D/C;

    double y_avg = std::accumulate(beg, y.end(), 0.0) / n_records;
    auto   exp_x = std::vector<double>(n_records, 0.0);

    for(int i = offset; i < n_areas; ++i) exp_x[i] = std::exp(-C*i);

    double exp_avg_x = std::accumulate(exp_x.begin(), exp_x.end(), 0.0)/n_records;

    double A = B * exp_avg_x + y_avg;

    return A;
  }
};

std::ostream& operator<<(std::ostream& o, trajectory t) {
  std::for_each(t.begin(), t.end(), [&o](PSpoint p) { o << p; });

  return o;
}

std::istream& operator>>(std::istream& i, trajectory& t) {
  PSpoint temp;
  while (i >> temp) {
    t.points.push_back(temp);
  }
  return i;
}

bool operator==(const trajectory& T, const trajectory& U) {
  int sizeT, sizeU;
  sizeT = T.size();
  sizeU = U.size();

  if (sizeT != sizeU) {
    return 0;
  }

  int equal_counter = 0;

  for (int i = 0; i < sizeT; ++i) {
    equal_counter += (T[i] == U[i]);
  }

  if (equal_counter == sizeT) {
    return 1;
  } else {
    return 0;
  }
}
