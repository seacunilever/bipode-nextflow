#!/usr/bin/env python

import argparse
import sys
import shutil
from pathlib import Path
import cmdstanpy
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def compile_stan_model(stan_file: Path) -> Path:
    """
    Compile a Stan model file.

    Args:
        stan_file: Path to the Stan model file

    Returns:
        Path to the compiled executable
    """
    logger.info(f"Compiling Stan model: {stan_file}")
    model = cmdstanpy.CmdStanModel(stan_file=stan_file)
    return Path(model.exe_file)

def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description='Compile a Stan model file')
    parser.add_argument('stan_file', type=Path, help='Path to the Stan model file')
    args = parser.parse_args()

    exe_file = compile_stan_model(args.stan_file)
    logger.info(f"Model compiled successfully: {exe_file}")
    sys.exit(0)

if __name__ == '__main__':
    main()
