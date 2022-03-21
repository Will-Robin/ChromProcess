from ChromProcess.Classes import SearchNode

values = [1,2,3,4,5]

root = SearchNode(1)

for v in values:
    root.insert_node(v)

def in_order_print(root):
    if root == None:
        return
    in_order_print(root.left)
    print(root.value)
    in_order_print(root.right)

def pre_order_print(root):
    if root == None:
        return        
    print(root.value)
    pre_order_print(root.left)
    pre_order_print(root.right)  

def find_value(root, value):
    if root == None:
        print("value not in tree")
    elif root.value == value:
        print(f"Found {value}")
    else:
        if root.value > value:
            find_value(root.left, value)
        else:
            find_value(root.right, value)

for v in values:
    find_value(root, v)

