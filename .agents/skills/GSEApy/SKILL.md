```markdown
# GSEApy Development Patterns

> Auto-generated skill from repository analysis

## Overview
This skill teaches the development patterns and conventions used in the GSEApy repository, a Rust-based codebase. It covers file naming, import/export styles, commit message conventions, and testing patterns. This guide is ideal for contributors aiming to maintain consistency and quality in their code contributions to GSEApy.

## Coding Conventions

### File Naming
- **Style:** camelCase
- **Example:**  
  ```plaintext
  geneSetAnalysis.rs
  enrichmentScore.rs
  ```

### Import Style
- **Style:** Relative imports
- **Example:**
  ```rust
  mod geneSetAnalysis;
  use crate::enrichmentScore::calculate_score;
  ```

### Export Style
- **Style:** Named exports
- **Example:**
  ```rust
  pub fn run_analysis() { ... }
  pub struct AnalysisResult { ... }
  ```

### Commit Messages
- **Type:** Conventional commits
- **Prefix:** `fix`
- **Average Length:** 67 characters
- **Example:**
  ```
  fix: correct calculation of enrichment scores in geneSetAnalysis
  ```

## Workflows

### Code Contribution
**Trigger:** When adding or updating code in the repository  
**Command:** `/contribute`

1. Create a new branch for your feature or bugfix.
2. Use camelCase for any new file names.
3. Use relative imports for modules.
4. Export functions and structs using named exports.
5. Write or update tests in files matching `*.test.*`.
6. Commit changes using the conventional commit format, prefixed with `fix` if applicable.
7. Open a pull request for review.

### Bug Fixing
**Trigger:** When fixing a bug in the codebase  
**Command:** `/bugfix`

1. Identify the bug and create a new branch.
2. Make the necessary code changes, following the coding conventions.
3. Update or add relevant tests in `*.test.*` files.
4. Commit your changes with a message like:
   ```
   fix: resolve issue with [describe bug briefly]
   ```
5. Push your branch and open a pull request.

## Testing Patterns

- **Framework:** Unknown (ensure to follow repository standards or ask maintainers)
- **Test File Pattern:** Files should match `*.test.*`
- **Example:**
  ```rust
  // enrichmentScore.test.rs
  #[test]
  fn test_calculate_score() {
      // test implementation
  }
  ```
- Place tests alongside or within the same module as the code being tested.

## Commands
| Command      | Purpose                                         |
|--------------|-------------------------------------------------|
| /contribute  | Start the code contribution workflow            |
| /bugfix      | Begin the bug fixing workflow                   |
```
