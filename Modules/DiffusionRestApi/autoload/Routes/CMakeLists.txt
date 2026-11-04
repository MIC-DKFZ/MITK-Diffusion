# MitkRESTAPI loads the libraries in its auto-load folder (lib/MitkRESTAPI/) when it is loaded: in
# every application with the REST server, and in none of the command line apps. This module
# registers the diffusion routes with the server service as soon as the service exists, so that
# the routes are there before any plug-in starts. A plug-in cannot do that: BlueBerry starts a
# plug-in that requires org.mitk.gui.qt.common only when one of its views is opened (an eager
# start would leave its views without their lifecycle notifications), and the Agent view looks
# the endpoints up before any diffusion view is open.
mitk_create_module(DiffusionRestApiRoutes
  DEPENDS PRIVATE MitkDiffusionRestApi
  AUTOLOAD_WITH MitkRESTAPI
)
