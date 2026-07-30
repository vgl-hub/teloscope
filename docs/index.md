<!--
This page reuses the repository README as the documentation home so there is a
single source of truth. The include-markdown plugin pulls in ../README.md at
build time and rewrites its relative links to work inside the docs site.
-->
{%
   include-markdown "../README.md"
   rewrite-relative-urls=true
%}
