setup_path=$(realpath ~/.inductiva/v0.4.4/completions/zsh/)
fpath=($setup_path $fpath)
autoload -Uz compinit
compinit
