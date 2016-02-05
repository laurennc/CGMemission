"""Runs tlac_web via bokeh server.

How to run it:
cd to_a_directotry_where_there_are_tlac_files
bokeh-server --script PATH_TO_HERE/tlac_web/bokeh_server.py

Then go to http://localhost:5006/tlac_web on your webbrowser.

To run it to be accessed from other computers do something like
bokeh-server --ip cursa.uio.no --port 12473 --script ../../../tlac_web.py

"""

from tlac_web import TlacWebApp
from bokeh.server.app import bokeh_app
from bokeh.server.utils.plugins import object_page


# The following code adds a "/tlac_web/" url to the bokeh-server.
@bokeh_app.route("/tlac_web/")
@object_page("tlac_web")
def make_object():
    app = TlacWebApp.create()    
    return app
    
