n = 100000
f = open('test_data/reallybig.in', 'w')
f.write('{"elements": [')
element = """
        {
            "pts": [
                [-20.0, 0.0],
                [20.0, 0.0]
            ],
            "bc_type": "traction",
            "bc": [
                [0.0, 0.0],
                [0.0, 0.0]
            ]
        }
"""
element_list = ','.join([element] * n)
f.write(element_list)
f.write(']}')
