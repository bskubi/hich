find . -name "*.py" | xargs awk '/^import[[:space:]]/ || /^from[[:space:]]/ { print FILENAME ": " $0 }'
