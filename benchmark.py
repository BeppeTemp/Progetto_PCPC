import matplotlib.pyplot as plt
result = dict()

with open('./strong.txt') as f1:
    for line in f1.readlines():
        key, *values = line.split(" ")
        values[-1] = values[-1].rstrip('\n')

        sum = 0
        for value in values:
            sum += float(value)

        result.update({key : {"value" : sum, "count" : 5,"efficeny":0,"speeedup":0}})

weak = dict()

with open('./weak.txt') as f1:
    for line in f1.readlines():
        key, *values = line.split(" ")
        values[-1] = values[-1].rstrip('\n')

        sum = 0
        for value in values:
            sum += float(value)

        size = f"{key * 100}x100"

        weak.update({key : {"value" : sum, "count" : 5,"efficeny":0,"problem_size":size}})

print("weak")
for key,value in weak.items():
    print(weak[key]["value"], weak[key]["count"])
    weak[key]["value"] /= weak[key]["count"]
    weak[key]["efficeny"] = weak["1"]["value"] / weak[key]["value"]
    print(f"{key}|{round(weak[key]['value'],2)}|{round(weak[key]['efficeny'],2)}|{weak[key]['problem_size']}")
print("strong")
for key,value in result.items():
    result[key]["value"] /= result[key]["count"]
    result[key]["efficeny"] = (result["1"]["value"] / (int(key) * result[key]["value"])) * 100
    result[key]["speedup"] = result["1"]["value"] / result[key]["value"]
    print(f"{key}|{round(result[key]['value'],2)}|{round(result[key]['efficeny'],2)}|{round(result[key]['speedup'],2)}")

plt.title("Strong Scalability")
plt.grid(True)
plt.xticks([x for x in range(1,17)])
#plt.yticks([x for x in range(1,17)])
plt.plot([x for x in range(1,17)],[result[key]["value"] for key in result],'o-')
plt.ylabel('Execution Time (sec)')
plt.xlabel('Number Processor')
plt.savefig("strong.png", dpi = 200)
plt.show()

plt.title("Strong Scalability Speedup")
plt.grid(True)
plt.xticks([x for x in range(1,17)])
plt.plot([x for x in range(1,17)],[x for x in range(1,17)])
plt.plot([x for x in range(1,17)],[result[key]["speedup"] for key in result],'o-')
plt.ylabel('speedup')
plt.xlabel('Number Processor')
plt.legend(['ideal speedup','real speedup'])
plt.savefig("speedup_strong.png", dpi = 200)
plt.show()


plt.title("Weak Scalability")
plt.grid(True)
plt.xticks([x for x in range(1,17)])
plt.plot([x for x in range(1,17)],[weak[key]["value"] for key in result],'o-')
plt.ylabel('Execution Time (sec)')
plt.xlabel('Number Processor')
plt.savefig("weak.png", dpi = 200)
plt.show()

plt.title("Weak Scalability Speedup")
plt.grid(True)
plt.xticks([x for x in range(1,17)])
plt.plot([x for x in range(1,17)],[1 for x in range(1,17)])
plt.plot([x for x in range(1,17)],[weak[key]["efficeny"] for key in result],'o-')
plt.ylabel('speedup')
plt.xlabel('Number Processor')
plt.legend(['ideal speedup','real speedup'])
plt.savefig("speedup_weak.png", dpi = 200)
plt.show()

