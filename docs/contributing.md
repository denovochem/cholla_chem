
1. Fork the repository
2. Clone your fork
3. Install dependencies
4. **Important**: Install pre-commit hooks: `pre-commit install`

## Development Process

1. Create a feature branch: `git checkout -b feature/your-feature`
2. Make your changes
3. Ensure all pre-commit hooks pass
4. Push your branch and create a Pull Request

## Code Quality

All code must pass:
- Ruff linting and formatting
- mypy type checking
- pytest test suite

Run these checks locally before pushing to ensure your PR will pass CI.

## Pre-commit Hooks

Pre-commit hooks are **required** for all contributors. They ensure:
- Consistent code style
- Type safety
- All tests pass before pushing
