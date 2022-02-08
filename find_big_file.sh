FILELINE=`git rev-list --all | xargs -rL1 git ls-tree -r --long | sort -uk3 | sort -rnk4 | head -10 | awk -F ' ' '{print $5}'| tr '\n' ' '`
echo $FILELINE
git filter-branch --force --tree-filter "rm -f $FILELINE" -- --all

