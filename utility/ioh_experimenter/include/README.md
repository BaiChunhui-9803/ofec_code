## Quickstart

### Installation

The following toolkits are needed for compiling IOHexperimenter:

* A `C++` compiler. The minimum compiler version is g++7 or equivalent, but we recommend g++9 or equivalent.
* [CMake](https://cmake.org), version 3.15 or higher

Please use the following commands to download, compile, and install this package:

```sh
> git clone --recursive https://github.com/IOHprofiler/IOHexperimenter.git
> cd IOHexperimenter
> mkdir build
> cd build
> cmake .. && make install
```

which installs all header files to `/usr/local/include/ioh` by default.
If you want to change this directory, please use the following flag `cmake -DCMAKE_INSTALL_PREFIX=your/path ..`

If you want to change build options, check the output of `cmake -L` or use `cmake-gui` or `ccmake`.

### Usage

To obtain a built-in problem, you could create a problem instance by passing the
*instance id* and the *search dimension* to the constructor of the intended problem class, e.g.,

```C++
#include "ioh.hpp"
auto om = ioh::problem::pbo::OneMax(1, 10); // PBO problem: instance 1, dim 10
auto sp = ioh::problem::bbob::Sphere(1, 5); // BBOB problem: instance 1, dim 5
```

The instance id is intended to generalize a certain problem by some transformations, where
it serves as the random seed for randomizing the transformations, e.g., affine
transforms for BBOB problems and scaling of objective values for PBO problems. Please see

* [PBO transformations](https://iohprofiler.github.io/IOHproblem/)
* [BBOB/COCO transformations](https://coco.gforge.inria.fr/downloads/download16.00/bbobdocfunctions.pdf)

Alternatively, we can also obtain a `std::shared_ptr` to a problem by using a problem registry of a given type. Using the registry, we can obtain problems by their name, which can be usefull in a dynamic setting:

```C++
const auto &problem_factory = ioh::problem::ProblemRegistry<ioh::problem::IntegerSingleObjective>::instance();
const std::shared_ptr<ioh::problem::IntegerSingleObjective> om = problem_factory.create("OneMax", 1, 10);
```

`Problem` objects are callable objects, and a search point can be evaluated on a problem by passing it to the call operator of the problem, for example:

```C++
std::vector<double> x0 = {1., 2., 3., 4., 5.};
double y0 = sp(x0);
```

For more advanced usage examples, please take a look at the example folder in the github repository. Some examples include:

* Using a [problem](./example/problem_example.hpp)
* Using a pre-defined [problem suite](./example/suite_example.hpp)
* Using the [logger](./example/logger_example.hpp) for storing benchmark data

### Adding new problems

We offer a very simple and convenient interface for integrating new benchmark problems/functions. First, you could define a new `test_problem` as you like. Note that the `<vector>` header is already imported in "ioh.hpp".

```C++
ioh::problem::SingleObjective test_problem(const std::vector<double> &)
{
    // the actual function body starts here
    // ...
    return 0.0;
}
```

Then, you only need to "wrap" this new function as follows:

```c++
ioh::problem::wrap_function<double, ioh::problem::SingleObjective>(
  &test_problem,
  "test_problem" // name for the new function
);
```

After wrapping, we could create this `test_problem` from the problem factory. Note that,
the instance id is ineffective in this approach since we haven't implemented any transformations for the wrapped problem.

```c++
auto &factory = ioh::problem::ProblemRegistry<ioh::problem::RealSingleObjective>::instance();
auto new_problem_f = factory.create(
  "test_problem",  // create by name
  1,               // instance id
  10               // number of search variables
);
```

Alternatively, one might wish to create the new problem by subclassing the abstract problem class
in IOHexperimenter, taking benefits of implementing more details, e.g., aforementioned transformations.
This can be done by inheriting the corresponding problem registration class, which is

* `ioh::problem::IntegerSingleObjective` for pseudo-Boolean problems, and
* `ioh::problem::RealSingleObjective` for continuous problems.

In the below example, we show how to do this for pseudo-Boolean problems.

```C++
class NewBooleanProblem final : public ioh::problem::RegisteredProblem<NewBooleanProblem, ioh::problem::IntegerSingleObjective>
{
protected:
    // [mandatory] The evaluate method is mandatory to implement
    double evaluate(const std::vector<int> &x) override
    {
        // the function body goes here
        return 0.0;
    }

    // [optional] If one wish to implement transformations on search variables
    std::vector<int> transform_variables(std::vector<int> x) override
    {
        return x;
    }

    // [optional] If one wish to implement transformations on objective values
    double transform_objectives(const double y) override
    { 
        return y;
    } 

public:
    /// [mandatory] This constructor always take `instance` as input even
    /// if it is ineffective by default. `instance` would be effective if and only if
    /// at least one of `transform_objectives` and `transform_objectives` is implemented
    NewBooleanProblem(const int instance, const int n_variables) :
        RegisteredProblem(
          ioh::problem::MetaData(
            1,                     // problem id, which will be overwritten when registering this class in all pseudo-Boolean problems
            instance,              // the instance id
            "NewBooleanProblem",   // problem name
            n_variables           // search dimensionality
            ) 
          )
    {
    }
};
```

Please check [this example](./example/problem_example.hpp) for adding continuous problems in this manner.
