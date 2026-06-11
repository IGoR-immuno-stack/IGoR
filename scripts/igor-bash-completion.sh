# bash completion for igor CLI
#
# source this file from ~/.bashrc or ~/.bash_profile:
#   source /path/to/igor-bash-completion.sh

_igor()
{
    local cur prev words cword

    # portable init: works with or without bash-completion package
    if declare -f _init_completion > /dev/null 2>&1; then
        _init_completion -n : || return
    else
        cword=$COMP_CWORD
        words=("${COMP_WORDS[@]}")
        cur="${COMP_WORDS[COMP_CWORD]}"
        prev="${COMP_WORDS[COMP_CWORD-1]}"
    fi

    local commands="datadir init config import-seqs align infer evaluate generate run"
    local global_opts="-h --help -v --version -w --workdir -b --batch -j --threads --stdout-file"
    local config_actions="get set list schema edit"
    local config_keys="
model.source
model.species
model.chain
model.parms
model.marginals
model.load_last_inferred
sequences.input
infer.iterations
infer.likelihood_threshold
infer.probability_ratio_threshold
infer.viterbi
infer.fix_error_rate
infer.only
infer.not
generate.seed
generate.fast
generate.error
generate.cdr3
generate.filename_prefix
generate.threads
pipeline.steps
"

    local command=""
    local i

    # find first non-option token = subcommand
    for (( i=1; i < cword; i++ )); do
        case "${words[i]}" in
            datadir|init|config|import-seqs|align|infer|evaluate|generate|run)
                command="${words[i]}"
                break
                ;;
            -w|--workdir|-b|--batch|-j|--threads|--stdout-file)
                ((i++))
                ;;
        esac
    done

    # complete values for global options
    case "$prev" in
        -w|--workdir|--stdout-file)
            compopt -o filenames 2>/dev/null
            COMPREPLY=( $(compgen -d -- "$cur") )
            return
            ;;
        --replay)
            compopt -o filenames 2>/dev/null
            COMPREPLY=( $(compgen -f -- "$cur") )
            return
            ;;
        -b|--batch)
            return
            ;;
        -j|--threads)
            COMPREPLY=( $(compgen -W "1 2 4 8 16" -- "$cur") )
            return
            ;;
        --gene)
            COMPREPLY=( $(compgen -W "V D J all" -- "$cur") )
            return
            ;;
    esac

    # no subcommand yet → complete subcommands or global options
    if [[ -z "$command" ]]; then
        if [[ "$cur" == -* ]]; then
            COMPREPLY=( $(compgen -W "$global_opts" -- "$cur") )
        else
            COMPREPLY=( $(compgen -W "$commands" -- "$cur") )
        fi
        return
    fi

    case "$command" in

        datadir|init|infer|evaluate)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$global_opts --help" -- "$cur") )
            fi
            return
            ;;

        import-seqs)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$global_opts --help" -- "$cur") )
            else
                compopt -o filenames 2>/dev/null
                COMPREPLY=( $(compgen -f -- "$cur") )
            fi
            return
            ;;

        align)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$global_opts --help --gene" -- "$cur") )
            fi
            return
            ;;

        generate)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$global_opts --help" -- "$cur") )
            else
                COMPREPLY=( $(compgen -W "1 10 100 1000" -- "$cur") )
            fi
            return
            ;;

        run)
            if [[ "$cur" == -* ]]; then
                COMPREPLY=( $(compgen -W "$global_opts --help --replay" -- "$cur") )
            else
                # propose manifests in .igor/runs/ of current workdir
                local workdir="."
                for (( i=1; i < cword; i++ )); do
                    if [[ "${words[i]}" == "-w" || "${words[i]}" == "--workdir" ]]; then
                        workdir="${words[i+1]}"
                        break
                    fi
                done
                local runs_dir="$workdir/.igor/runs"
                if [[ -d "$runs_dir" ]]; then
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "$runs_dir/$cur") )
                else
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "$cur") )
                fi
            fi
            return
            ;;

        config)
            local action=""
            for (( i=1; i < cword; i++ )); do
                case "${words[i]}" in
                    get|set|list|schema|edit)
                        action="${words[i]}"
                        break
                        ;;
                esac
            done

            if [[ -z "$action" ]]; then
                if [[ "$cur" == -* ]]; then
                    COMPREPLY=( $(compgen -W "$global_opts --help" -- "$cur") )
                else
                    COMPREPLY=( $(compgen -W "$config_actions" -- "$cur") )
                fi
                return
            fi

            case "$action" in
                get)
                    COMPREPLY=( $(compgen -W "$config_keys" -- "$cur") )
                    return
                    ;;
                set)
                    local seen_key=""
                    for (( i=1; i < cword; i++ )); do
                        if [[ "${words[i]}" == "set" && $((i+1)) -lt cword ]]; then
                            seen_key="${words[i+1]}"
                            break
                        fi
                    done

                    if [[ "$prev" == "set" ]]; then
                        COMPREPLY=( $(compgen -W "$config_keys" -- "$cur") )
                    else
                        case "$seen_key" in
                            model.source)
                                COMPREPLY=( $(compgen -W "custom builtin last_inferred" -- "$cur") )
                                ;;
                            model.load_last_inferred|infer.viterbi|infer.fix_error_rate|\
                            generate.fast|generate.error|generate.cdr3)
                                COMPREPLY=( $(compgen -W "true false" -- "$cur") )
                                ;;
                            model.parms|model.marginals|sequences.input)
                                compopt -o filenames 2>/dev/null
                                COMPREPLY=( $(compgen -f -- "$cur") )
                                ;;
                            pipeline.steps)
                                COMPREPLY=( $(compgen -W "import-seqs align infer evaluate generate" -- "$cur") )
                                ;;
                            model.species)
                                COMPREPLY=( $(compgen -W "human mouse" -- "$cur") )
                                ;;
                            model.chain)
                                COMPREPLY=( $(compgen -W "alpha beta delta gamma" -- "$cur") )
                                ;;
                            *)
                                COMPREPLY=()
                                ;;
                        esac
                    fi
                    return
                    ;;
                list|schema|edit)
                    if [[ "$cur" == -* ]]; then
                        COMPREPLY=( $(compgen -W "$global_opts --help" -- "$cur") )
                    fi
                    return
                    ;;
            esac
            ;;
    esac
}

complete -F _igor igor
