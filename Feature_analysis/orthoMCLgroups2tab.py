#!/usr/bin/python

import sys

class Cluster:
    def __init__(self, name):
        self.name = name
        self.members = {}

    def add(self, strain, member):
        self.members[strain] = member

all_proteins = {}
clusters = []
strains = {}

for ln in open(sys.argv[1]):
    if ln.startswith('>'):
        protein = ln[1:].rstrip()
        all_proteins[protein] = False

n = 0
for ln in open(sys.argv[2]):
    prefix, proteins = ln.rstrip().split(':')
    proteins = proteins.split()

    c = Cluster(prefix)

    for p in proteins:
        strain, id = p.split("|")
        c.add(strain, id)
        strains[strain] = strain

        all_proteins[p] = True

    clusters.append(c)

# add the singletones
i = 1
for p, already_counted in all_proteins.iteritems():
    if already_counted == False:
        c = Cluster("single%d" % (i))
        strain, id = p.split("|")
        c.add(strain, id)
        i += 1
        clusters.append(c)


first_cluster = clusters[0]
print " ".join(['"%s"' % (c.name) for c in clusters])

i = 0
for b in sorted(strains):
    print '"%s"' % (b),

    for c in clusters:
        print " ",
        if b in c.members:
            print "1",
        else:
            print "0",
    print

    i += 1
