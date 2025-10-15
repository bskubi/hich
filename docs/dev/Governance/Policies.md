# Code changes
1. Create a GitHub issue using the "code change" template.
2. Create a git branch named `I[num]-[abbr_desc]`
3. Explicitly define
	1. Scope
	2. Definition of done
# Questions
* Should we have separate "dev" and "main" branches?
* What if highly related work outside pre-defined scope is discovered while working on an issue?

# Scope creep
+ If new related issues are identified during work on an issue
1. Create new GitHub issue immediately
2. Implement out-of-scope changes as a separate issue+branch.
3. Crosslink in new GitHub issue with issue where it was discovered.

# Storing information about development process

|                             |                                                                                                                                                       |                                                                                                                                                                                              |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Location                    | Purpose & Audience                                                                                                                                    | Content from Your Log                                                                                                                                                                        |
| 📝 GitHub Issue             | Public Record of the Problem & Solution. For users and your future self to understand what the problem was and that it has been resolved.             | • **Initial Post:** The problem description.  <br>• **Final Comment:** A concise summary of the fix. Link to the commit(s) and dev log.                                                      |
| 💻 Git Commit Message       | Code-Level Explanation. For a developer reading the code history (git log) to understand why a specific change was made. It should be self-contained. | • **Subject Line:** A short imperative summary.<br>• **Body:** A paragraph explaining the problem and the implemented solution. Format: cause of problem, change made to implement solution. |
| 📓 Developer Log (Your Doc) | The Full Story. Your personal/internal "lab notebook" for the code. Useful to reconstruct your thought process, including dead ends and future ideas. | The entire, unabridged log.                                                                                                                                                                  |
