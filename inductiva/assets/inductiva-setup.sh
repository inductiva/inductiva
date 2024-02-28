script_dir=$(realpath "$(dirname "$0")")
setup_path=$(realpath "$script_dir/completions/zsh/")
fpath=($setup_path $fpath)
autoload -Uz compinit
compinit
