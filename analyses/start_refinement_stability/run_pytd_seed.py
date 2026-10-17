#!/usr/bin/env python3
"""Run PyTransDecoder (longorfs + predict) with the start-PWM training RNG seeded by --seed.

PyTransDecoder seeds start-PWM training with random.Random(1). This driver substitutes
another seed without modifying the package, so the sensitivity of start refinement
to that seed can be measured.

usage: run_pytd_seed.py --seed N -- <pytransdecoder cli args>   (run once for longorfs, once for predict)
"""
import random
import sys

argv = sys.argv[1:]
seed = int(argv[argv.index("--seed") + 1])
cli_args = argv[argv.index("--") + 1:]

import pytransdecoder.predict as predict  # noqa: E402


class _SeededRandom(random.Random):
    def __init__(self, _ignored=None):
        super().__init__(seed)


class _RandomModule:
    def __getattr__(self, name):
        return _SeededRandom if name == "Random" else getattr(random, name)


predict.random = _RandomModule()

from pytransdecoder import cli  # noqa: E402

sys.argv = ["pytransdecoder"] + cli_args
sys.exit(cli.main())
