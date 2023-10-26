import re

if __name__ == "__main__":
    with open("data/corona.txt", "r") as f:
        s = re.sub(r"\s+", '', f.read())
        print(type(s), s)
