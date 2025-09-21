# AI Coding Agent Instructions

## Project Overview
This is a numerical analysis homework project implementing Newton's method for root finding in C++. The codebase contains a single file `newton.cpp` with an abstract `Function` class and an incomplete `NewtonSolver` class.

## Key Architecture

### Core Components
- **`Function` class**: Abstract base class defining mathematical functions
  - `operator()`: Evaluates f(x) 
  - `d()`: Evaluates the derivative f'(x)
- **`NewtonSolver` class**: Implements Newton's iterative root-finding algorithm
  - Takes a `Function` reference in constructor
  - `solve(double x0)`: Main method to implement - finds root starting from x0

### Implementation Requirements
The `solve()` method must implement Newton's iteration formula: `x_{n+1} = x_n - f(x_n)/f'(x_n)`

## Development Patterns

### Code Style
- Uses object-oriented design with inheritance
- Constructor takes function by reference (not pointer)
- Comments are in Chinese (数值分析 context)
- Method implementations should be concise and focused

### Testing Approach
- Test code is provided in comments but should NOT be included in submissions
- The example `Func1` class implements `f(x) = x - tan(x)` with derivative `f'(x) = 1 - sec²(x)`
- Test starting point is x0 = 4.5

### Submission Guidelines
- **Critical**: Remove all test code (main function, Func1 class) before submission
- Only submit the core `Function` and `NewtonSolver` class implementations
- The solve method should handle convergence and return the final root approximation

## Common Implementation Considerations
- Add convergence criteria (tolerance and max iterations)
- Handle division by zero when f'(x) ≈ 0
- Consider numerical stability for the iteration
- Return type is `double` - ensure proper precision

## Build Instructions
Standard C++ compilation - no special build system required:
```bash
g++ -o newton newton.cpp
```