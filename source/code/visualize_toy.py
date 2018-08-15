import sys

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph 

username='aritz@reed.edu'
password='platypus'



def visualize_graph(nxG,pos,neg,predictions,title):

    gsG = GSGraph()
    gsG.set_name(title)
    gsG.set_tags(['CREU'])

    NEG_COLOR = rgb2hex(1,.5,0,normed=True)
    POS_COLOR = rgb2hex(0,.5,1,normed=True)

    for node in nxG.nodes:
        print(node,predictions[node])
        if 'prime' in node:
            gsG.add_node(node,label=node.replace('_prime','\''),popup='Prediction: %.4f' % (predictions[node]))
            size=30
        else:
            gsG.add_node(node,label=node.replace('_',' L'),popup='Prediction: %.4f' % (predictions[node]))
            size=20

        if node in pos:
            color = POS_COLOR
            shape='star'
        elif node in neg:
            color = NEG_COLOR
            shape='star'
        else: # unlabeled
            pred = predictions[node]
            color = rgb2hex(1-pred,0.5,pred,normed=True)
            shape='rectangle'
        gsG.add_node_style(node, color=color, shape=shape, height=size, width=size*2)

    for n1,n2 in nxG.edges:
        gsG.add_edge(n1,n2,popup='Weight: %.4f' % (nxG[n1][n2]['weight']))
        gsG.add_edge_style(n1,n2,width=3*nxG[n1][n2]['weight'])

    gid = post(gsG,username,password,group='CREU')
    print('Graph posted with id',gid)

    return


def post(G,username,password,group=None):
  '''
  Post a graph to graphspace.  If this graph is new, share it with
  the group.
  Inputs: Graph (GSGraph Object)
  Outputs: the graph ID (int)
  '''
  print('Posting graph...')
  
  # connect to GraphSpace with the username and password.
  gs = GraphSpace(username,password)
  try:
    # try updating the graph. If the graph does not exist, 
    # this will throw an error.
    graph = gs.update_graph(G)
  except:
    # catch the error and try posting a new graph.
    graph = gs.post_graph(G)
    if group:
        # share this graph with the group.
        gs.share_graph(graph_id=graph.id,group_name=group)
  return graph.id

#http://www.psychocodes.in/rgb-to-hex-conversion-and-hex-to-rgb-conversion-in-python.html
def rgb2hex(r,g,b,normed=False):
    if normed:
        hexv = "#{:02x}{:02x}{:02x}".format(int(r*255),int(g*255),int(b*255))
    else:
        hexv = "#{:02x}{:02x}{:02x}".format(r,g,b)
    return hexv

def hex2rgb(hexcode):
    rgb = tuple(map(ord,hexcode[1:].decode('hex')))
    return rgb


if __name__ == '__main__':
    main(sys.argv)



            




