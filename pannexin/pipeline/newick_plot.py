#!/usr/bin/env python3
"""Minimal Newick → rectangular phylogram (no Biopython required)."""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class Node:
    name: str | None = None
    length: float = 0.0
    support: float | None = None
    children: list[Node] = field(default_factory=list)
    # layout
    x: float = 0.0
    y: float = 0.0


def _tokenize(newick: str) -> list[str]:
    s = newick.strip()
    if s.endswith(";"):
        s = s[:-1]
    return re.findall(r"[(),]|[^(),]+", s)


def parse_newick(newick: str) -> Node:
    tokens = _tokenize(newick)
    pos = 0

    def peek() -> str | None:
        return tokens[pos] if pos < len(tokens) else None

    def take() -> str:
        nonlocal pos
        tok = tokens[pos]
        pos += 1
        return tok

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if peek() == "(":
            take()
            node.children.append(parse_node())
            while peek() == ",":
                take()
                node.children.append(parse_node())
            if peek() != ")":
                raise ValueError("expected ')' in Newick")
            take()
            # optional internal label / bootstrap
            if peek() and peek() not in {",", ")", ":", "("}:
                label = take().strip()
                if ":" in label:
                    lab, _, rest = label.partition(":")
                    label, length_s = lab, rest
                    try:
                        node.length = float(length_s)
                    except ValueError:
                        pass
                if label:
                    try:
                        node.support = float(label)
                    except ValueError:
                        node.name = label
        # tip or length after internal
        if peek() and peek() not in {",", ")", "("}:
            label = take().strip()
            if ":" in label:
                name, _, length_s = label.partition(":")
                if name and node.name is None and node.support is None:
                    # tip name
                    if not node.children:
                        node.name = name
                    else:
                        try:
                            node.support = float(name)
                        except ValueError:
                            node.name = name
                try:
                    node.length = float(length_s)
                except ValueError:
                    pass
            else:
                if not node.children:
                    node.name = label
                else:
                    try:
                        node.support = float(label)
                    except ValueError:
                        node.name = label
        elif peek() == ":":
            take()
            if peek():
                try:
                    node.length = float(take())
                except ValueError:
                    pass
        return node

    root = parse_node()
    return root


def ladderize(node: Node) -> int:
    """Sort children by subtree size; return tip count."""
    if not node.children:
        return 1
    sizes = []
    for ch in node.children:
        sizes.append((ladderize(ch), ch))
    sizes.sort(key=lambda t: t[0])
    node.children = [ch for _, ch in sizes]
    return sum(s for s, _ in sizes)


def layout(node: Node) -> None:
    """Assign x = distance from root, y = tip index."""
    tips: list[Node] = []

    def collect(n: Node) -> None:
        if not n.children:
            tips.append(n)
        else:
            for ch in n.children:
                collect(ch)

    collect(node)

    def set_y(n: Node) -> float:
        if not n.children:
            return n.y
        ys = [set_y(ch) for ch in n.children]
        n.y = sum(ys) / len(ys)
        return n.y

    for i, tip in enumerate(tips):
        tip.y = float(i)
    set_y(node)

    def set_x(n: Node, dist: float) -> None:
        n.x = dist
        for ch in n.children:
            set_x(ch, dist + max(ch.length, 0.0))

    set_x(node, 0.0)


def iter_edges(node: Node):
    for ch in node.children:
        yield node, ch
        yield from iter_edges(ch)


def iter_tips(node: Node):
    if not node.children:
        yield node
    else:
        for ch in node.children:
            yield from iter_tips(ch)
