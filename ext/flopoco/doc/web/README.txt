To produce the tgz of a version:
first edit VERSION,
then in bash:
git archive --prefix=flopoco-$(cat VERSION)/ -o flopoco-$(cat VERSION).tgz HEAD
Place it on the forge;
Then edit  flopoco_installation_manual.html using the path created by the forge.
commit and push, and scp it as below.
Finally do a 
git tag -a flopoco-$(cat VERSION) -m"Release!"

git push --tags


To update the web site from an admin account: simply commit and push in the web directory of the master branch.
gitlab will transfer the pages to the web site.
Instructions to manage this:
- The script to run is in .gitlab-ci.yml
- Florent has a worker on https://ci.inria.fr/
- for the rest see https://inria-ci.gitlabpages.inria.fr/doc/page/web_portal_tutorial/

The operators.html page is generated from flopoco itself:
./flopoco BuildHTMLDoc
(due to lazy hardcoding of relative paths, it must be run from the base flopoco directory, i.e. not from the build/ directory. Same for autotest) 
THIS IS CURRENTLY NOT DONE IN CI because we don't re-compile flopoco in CI yet.


The html bibliography pages are generated using bibtex2html
1/ edit doc/web/bib/flopoco.bib or doc/web/bib/flopoco-users.bib
2/ commit and push.

This should lauch the commands previously managed by hand:
BEGIN HISTORY now in .gitlab-ci.yml
In the doc/web/bib directory, run:
bibtex2html -t "Publications about FloPoCo" --header "<p><em>If some of your works belong there, please drop a mail to F. de Dinechin with the corresponding bibtex entries</em></p><hr>" -d -r -revkeys flopoco.bib
bibtex2html -t "Publications using FloPoCo" --header "<p><em>If some of your works belong there, please drop a mail to F. de Dinechin with the corresponding bibtex entries</em></p><hr>" -d -r -revkeys flopoco-users.bib
END HISTORY

3
