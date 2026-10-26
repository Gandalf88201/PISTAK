import sys

from pistak.cli import main

# Il controllo serve a multiprocessing (elaborazione in parallelo): i processi figli reimportano questo modulo.
if __name__ == "__main__":
    sys.exit(main())
