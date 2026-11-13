// Posts build/coverage/comment.md (scripts/tests/coverage_markdown.py) on the pull request, as
// one comment that every run edits in place. ci.yml calls it through actions/github-script,
// which passes the authenticated client and the event's context.
const fs = require('fs');

// Opens the comment, so that the next run finds it among the pull request's comments.
const MARKER = '<!-- igor-coverage -->';

module.exports = async ({ github, context }) => {
    const run = `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId}`;
    const body = [
        MARKER,
        fs.readFileSync('build/coverage/comment.md', 'utf8'),
        `The HTML report is the coverage-report artifact of [this run](${run}).`,
    ].join('\n');

    const comments = await github.paginate(github.rest.issues.listComments, {
        ...context.repo,
        issue_number : context.issue.number,
        per_page : 100,
    });
    const previous = comments.find((c) => c.user?.type === 'Bot' && c.body?.startsWith(MARKER));
    if (previous) {
        await github.rest.issues.updateComment({...context.repo, comment_id : previous.id, body });
    } else {
        await github.rest.issues.createComment({...context.repo, issue_number : context.issue.number, body });
    }
};
