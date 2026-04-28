import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1], label="Test")
leg = fig.legend(bbox_to_anchor=(0.5, 0.995))
print(fig.legends[0].get_bbox_to_anchor())
fig.legends[0].set_bbox_to_anchor((0.5, 0.95))
print(fig.legends[0].get_bbox_to_anchor())
