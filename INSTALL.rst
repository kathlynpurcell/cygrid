Installation procedure for Cygrid
=================================

Requirements
------------

Cygrid depends on the numpy and astropy Python libraries. If you want to install from source, you'll need Cython as well.

Quick installation with Pip
---------------------------

The quickest way to install Cygrid is with `pip <http://www.pip-installer.org>`_
(>= 1.4.2), which automatically fetches the latest version of Cygrid and any
missing dependencies::

    pip install --user cygrid

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade cygrid

Almost-as-quick installation from official source release
---------------------------------------------------------

Cygrid is also available in the
`Python Package Index (PyPI) <https://pypi.python.org/pypi/cygrid>`_. You can
download it with::

    curl -O https://pypi.python.org/packages/source/h/cygrid/cygrid-0.1.0.tar.gz

and build it with::

    tar -xzf cygrid-0.1.0.tar.gz
    pushd cygrid-0.1.0
    python setup.py install --user
    popd
