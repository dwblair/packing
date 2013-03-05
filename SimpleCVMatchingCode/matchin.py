from SimpleCV import Image

for i in range(10):
    template = Image('wave1.jpg')
    img = Image('wave2.jpg')
    match = img.drawKeypointMatches(template,490.00,0.05)
    match.save("match.jpg")
    print i
