Two primary branches: `main` and `dev`

Minor refactorings and fixes (individually < 1 day):
- Can work on multiple small sequential changes as atomic commits on `dev` branch.
Major, time-consuming, complex features and restructuring (> 1 day):
+ Branch from dev for each major feature.
Merge `dev` into `main` with a batch commit message once per week. Optionally, tag it as a new release.