n = 10000
f = open('data/reallybig.in', 'w')
f.write('{"elements": [')
element = """{"pts": [ [-20.0, 0.0],
                [20.0, 0.0]
            ],
            "bc_type": "traction",
            "bc": [
                [0.0, 0.0],
                [0.0, 0.0]
            ],
            "refine": 0
        }
""".replace(' ', '')
element_list = ','.join([element] * n)
f.write(element_list)
f.write(']}')
