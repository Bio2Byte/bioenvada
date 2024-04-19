
#!/usr/bin/python3
from ete3 import Tree
import sys


tree = sys.argv[1]
plot = sys.argv[2]

form = 0
if plot != 'true':
    form = int(plot)


t=Tree(tree, format=form)
image_name = tree.split('.')[0]+"_tree.png"

t.render(image_name, dpi=300)