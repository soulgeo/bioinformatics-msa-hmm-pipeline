from tabulate import tabulate

def print_float_table(data):
    first_key = next(iter(data))
    headers = ["ID"] + list(data[first_key].keys())
    table_rows = [[k] + list(v.values()) for k, v in data.items()]
    
    print(tabulate(table_rows, headers=headers, tablefmt="plain", floatfmt=".2f"))
