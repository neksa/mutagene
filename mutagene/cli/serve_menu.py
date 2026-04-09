"""CLI command for starting the local web server."""

import logging

logger = logging.getLogger(__name__)


class ServeMenu:
    def __init__(self, parser):
        parser.add_argument(
            "--host",
            default="127.0.0.1",
            help="Host to bind to (default: 127.0.0.1)",
        )
        parser.add_argument(
            "--port",
            type=int,
            default=5000,
            help="Port to bind to (default: 5000)",
        )
        parser.add_argument(
            "--no-browser",
            action="store_true",
            help="Do not automatically open browser",
        )
        parser.add_argument(
            "--debug",
            action="store_true",
            help="Enable debug mode",
        )

    def callback(self, args):
        import signal
        import sys

        try:
            from mutagene.webapp.server import start_server
        except ImportError:
            print(
                "Web application dependencies are not installed.\n"
                "Install them with: pip install mutagene[web]"
            )
            sys.exit(1)

        # Reset SIGINT handler to default for Flask-SocketIO
        # The main app's signal handler interferes with proper shutdown
        signal.signal(signal.SIGINT, signal.default_int_handler)

        logger.info(f"Starting web server at http://{args.host}:{args.port}")
        try:
            start_server(
                host=args.host,
                port=args.port,
                debug=args.debug,
                open_browser=not args.no_browser,
            )
        except KeyboardInterrupt:
            print("\n\nShutting down server...")
            sys.exit(0)
