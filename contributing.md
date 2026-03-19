# Contributing to omnomnomics

Thank you for contributing to **omnomnomics**. This project is developed collaboratively within the Prange Lab, with an emphasis on reproducibility, maintainability, and clear workflow behavior on HPC systems.

The guidelines below are intentionally lightweight.

---

## General workflow

1. **Create a feature branch** from `main`
   ```bash
   git checkout -b my-feature
   ```

2. **Make your changes**
   - Keep commits focused and readable
   - Prefer small, logically grouped commits over large ones

3. **Open a pull request (PR)** into `main`
   - Describe *what* you changed and *why*
   - Reference relevant issues or discussions if applicable

4. **Wait for review**
   - At least one approval is required before merging
   - Address comments and requested changes

Direct pushes to `main` are intentionally restricted.

---

## Coding guidelines

- Prefer clear, explicit code over clever abstractions
- Keep module responsibilities narrow
- Preserve the modular workflow design of the pipeline
- Keep comments factual and tied to the code
- Avoid overly long comments
- Avoid unnecessary conversational or self-referential comments in code

Where applicable:
- Keep CLI behavior stable unless a change is intentional
- Keep workflow file contracts and stage boundaries easy to reason about
- Prefer reproducible, inspectable intermediate outputs over opaque shortcuts

---

## Configuration and workflow behavior

- Changes to defaults or interfaces should be clearly documented in the PR
- Prefer extending structured configuration over introducing ad-hoc flags
- Keep workflow, site, and per-run configuration concerns separate
- If user-facing stage behavior changes, update the relevant documentation

---

## Testing

- New features should include basic validation where feasible
- Tests and checks should be lightweight enough to run in development environments
- If full workflow testing is too expensive, validate the smallest meaningful unit
- When changing workflow contracts, verify expected outputs and input assumptions explicitly

CI coverage may remain limited for some time, so local verification and focused tests are important.

---

## Documentation

- Update `README.md` if user-facing behavior changes
- Keep documentation concise and factual
- Do not describe unfinished or experimental features as stable defaults
- Keep `DESIGNGOALS.md` aligned with major design decisions during active redesign work

---

## Scope and stability

Parts of the pipeline are still being actively redesigned.

If you are contributing to active redesign areas:
- expect interfaces and file layouts to change
- keep backward compatibility only when it is genuinely useful
- document assumptions and constraints clearly in the PR

---

## Questions

If you are unsure about design decisions, scope, or implementation details:
- open a draft PR
- or discuss informally before investing significant effort
