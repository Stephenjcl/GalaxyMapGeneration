#Viability test
from garnets import generate_stellar_system, random_star
from enviroment import BreathabilityPhrase

breathableair = False
moon = False
attempts = 0

# while breathableair is False or attempts <= 35000:
#     try:
#         stellar_system = generate_stellar_system(random_star())
#         for i in range(0, len(stellar_system.planets)):
#             if stellar_system.planets[i].breath == BreathabilityPhrase.BREATHABLE:
#                 breathableair = True
#             for n in range(0, len(stellar_system.planets[i].moons)):
#                 if stellar_system.planets[i].moons[n].breath == BreathabilityPhrase.BREATHABLE:
#                     breathableair = True
#                     moon = True
#         attempts += 1
#
#         print('\r', attempts, "attempts made.", end='', flush=True)
#
#     except:
#         pass
#
# if moon == True:
#     type = "moon"
# else:
#     type = "planet"
# print('\r', "Habitable", type, "found after", attempts, "attempts. This gives a breathable atmosphere probability of:", 1/attempts*100, "%", end = '', flush=True)

#Count breathable atmospheres

# f = open('Production Run.csv', 'w')
#
# f.close()

with open('Production Run.csv', 'r') as content_file:
    content = content_file.read()
    print(content.count("BreathabilityPhrase.BREATHABLE"))
    habitable = content.count("BreathabilityPhrase.BREATHABLE")

print("There are", habitable, "habitable worlds in this galaxy. This means there is oxygen present, and no excess poisonous gases.", '\n',
"This means there is a", habitable/35000*100, "% chance of breathing the air on a given world.")
