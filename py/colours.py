import pandas as pd
from pastml.visualisation.colour_generator import get_enough_colours


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', nargs='+', type=str)
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()

    states = set()
    for f in params.input:
        states |= set(pd.read_csv(f, sep='\t')['consensus'].unique())
    states = sorted([_ for _ in states if 'sensitive' != _])
    one_mutation_states = sorted({_.split('+')[0] for _ in states})
    colours = get_enough_colours(len(one_mutation_states))
    state2colour = dict(zip(one_mutation_states, colours))

    with open(params.output, 'w+') as f:
        f.write('state\tcolour\n')
        f.write('sensitive\t#909090\n')
        for state in states:
            one_mut_state = state.split('+')[0]
            f.write('{}\t{}\n'.format(state, state2colour[one_mut_state]))
