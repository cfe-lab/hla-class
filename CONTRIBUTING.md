This is a Django app deployed on the server with some other Django apps.

## Deploying a Release
You need to connect to the server from inside the lab network. You can find the
internal IP address in the `servers.md` file in the `dev-docs` repository. If
you don't already have an account on the server, ask the lab director for
access.

The source code is deployed with a Git clone at
`/alldata/hla_class`. To deploy an update, change to that
directory, and then run `sudo git pull`. 

Note that this tool is not part of the Django framework, as it is written in
Ruby and PHP. There is an alias on the server which redirects from the 
Django web address to the above folder containing this tool.

Then restart the Apache server with
`sudo systemctl restart httpd`