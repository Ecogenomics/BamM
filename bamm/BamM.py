#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamM.py                                                         #
#                                                                             #
#    Description!!                                                            #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TemplateClass():
    """Utilities wrapper"""
    def __init__(self): pass

    def sayHi(self):
        print('write some "REAL" code you bum!')

    def demoStuff(self):

        """
        # parse a file
        try:
            with open(filename, "r") as fh:
                for line in fh:
                    print line
        except:
            print "Error opening file:", filename, exc_info()[0]
            raise
        """

        """
        fig = plt.figure()

        #-----
        # make a 3d plot
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(points[:,0],
                   points[:,1],
                   points[:,2],
                   #edgecolors='none',
                   #c=colors,
                   #s=2,
                   #marker='.'
                   )

        #-----
        # make a 2d plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(points[:,0],
                points[:,1],
                '*g')

        #-----
        # show figure
        plt.show()
        # or save figure
        plt.savefig(filename,dpi=300,format='png')

        #-----
        # clean up!
        plt.close(fig)
        del fig
        """

        return 0
