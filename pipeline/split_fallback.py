"""Split the browser fallback output (saved tool-result JSON) into yr_all.csv and bw_all.csv."""
import json, sys
d = json.load(open(sys.argv[1])); t = d[0]['text'] if isinstance(d, list) else d
if t.startswith('"'): t = json.loads(t[:t.rindex('"')+1])
t = t.split('\n#')[0]
yr, bw = t.split('@@yr\n')[1].split('\n@@bw\n')
open('yr_all.csv', 'w').write(yr + '\n'); open('bw_all.csv', 'w').write(bw + '\n')
print('yr lines', yr.count('\n'), 'bw lines', bw.count('\n'))
