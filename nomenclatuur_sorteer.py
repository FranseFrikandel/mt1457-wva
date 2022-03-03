import re

nc_lines = []
# String van characters om uit te sorteren voor bepalen van alfabetische volgorde
substr = "\ |\?|\.|\!|\/|\;|\:|\{|\}|\_|\$|\&(.*)"

with open("nomenclatuur_ongesorteerd.txt") as nc_unsorted:
    for line in nc_unsorted.readlines():
        nc_line_sort = re.sub(substr, "", line)
        print(nc_line_sort)
        nc_lines.append([line, nc_line_sort])

nc_lines.sort(key=lambda li : li[1])

with open("nomenclatuur_gesorteerd.txt", "w") as nc_sorted:
    for line in nc_lines:
        nc_sorted.write(line[0])
